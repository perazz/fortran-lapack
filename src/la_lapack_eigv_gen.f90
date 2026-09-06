!> Nonsymmetric eigenvalue, Schur and generalized Schur drivers
module la_lapack_eigv_gen
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_eigv_comp
     use la_lapack_eigv_comp2
     use la_lapack_eigv_gen2
     use la_lapack_eigv_gen3
     use la_lapack_eigv_gen_hess
     use la_lapack_givens_jacobi_rot
     use la_lapack_orthogonal_factors_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgges
     public :: la_sggesx
     public :: la_sggev
     public :: la_sggevx
     public :: la_sgees
     public :: la_sgeesx
     public :: la_sgeev
     public :: la_sgeevx
     public :: la_sgges3
     public :: la_sggev3
     public :: la_dgges
     public :: la_dggesx
     public :: la_dggev
     public :: la_dggevx
     public :: la_dgees
     public :: la_dgeesx
     public :: la_dgeev
     public :: la_dgeevx
     public :: la_dgges3
     public :: la_dggev3
     public :: la_qgges
     public :: la_qggesx
     public :: la_qggev
     public :: la_qggevx
     public :: la_qgees
     public :: la_qgeesx
     public :: la_qgeev
     public :: la_qgeevx
     public :: la_qgges3
     public :: la_qggev3
     public :: la_cgges
     public :: la_cggesx
     public :: la_cggev
     public :: la_cggevx
     public :: la_cgees
     public :: la_cgeesx
     public :: la_cgeev
     public :: la_cgeevx
     public :: la_cgges3
     public :: la_cggev3
     public :: la_zgges
     public :: la_zggesx
     public :: la_zggev
     public :: la_zggevx
     public :: la_zgees
     public :: la_zgeesx
     public :: la_zgeev
     public :: la_zgeevx
     public :: la_zgges3
     public :: la_zggev3
     public :: la_wgges
     public :: la_wggesx
     public :: la_wggev
     public :: la_wggevx
     public :: la_wgees
     public :: la_wgeesx
     public :: la_wgeev
     public :: la_wgeevx
     public :: la_wgges3
     public :: la_wggev3

     contains

     !> SGGES: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> SGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_sgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_s) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,maxwrk,minwrk
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'SGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'SORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'SORGQR',' ',n,1,n, &
                              -1))
                 end if
              else
                 minwrk = 1
                 maxwrk = 1
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -19
           end if
           if (info /= 0) then
              call la_xerbla('SGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           safmin = la_slamch('S')
           safmax = one/safmin
           call la_slabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n space for storing balancing factors)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_sggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_slaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_sorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_slaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_sgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_shgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 40
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: need 4*n+16 )
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_stgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_sggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_sggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_slascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           40 continue
           work(1) = maxwrk
           return
     end subroutine la_sgges
     !> DGGES: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> DGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_dgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_d) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,maxwrk,minwrk
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'DGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'DORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'DORGQR',' ',n,1,n, &
                              -1))
                 end if
              else
                 minwrk = 1
                 maxwrk = 1
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -19
           end if
           if (info /= 0) then
              call la_xerbla('DGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           safmin = la_dlamch('S')
           safmax = one/safmin
           call la_dlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n space for storing balancing factors)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_dggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_dlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_dorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_dlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_dgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_dhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 50
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: need 4*n+16 )
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_dtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_dggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_dggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_dlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           50 continue
           work(1) = maxwrk
           return
     end subroutine la_dgges
     !> QGGES: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> QGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_qgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_q) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,maxwrk,minwrk
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'QGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'QORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'QORGQR',' ',n,1,n, &
                              -1))
                 end if
              else
                 minwrk = 1
                 maxwrk = 1
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -19
           end if
           if (info /= 0) then
              call la_xerbla('QGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           safmin = la_qlamch('S')
           safmax = one/safmin
           call la_qlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n space for storing balancing factors)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_qggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_qlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_qorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_qlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_qgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_qhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 50
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: need 4*n+16 )
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_qtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_qggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_qggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_qlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           50 continue
           work(1) = maxwrk
           return
     end subroutine la_qgges

     !> SGGESX: computes for a pair of N-by-N real nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the real Schur form (S,T), and,
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_sggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim, &
     alphar,alphai,beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,iwork,liwork, &
               bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),rconde(2),rcondv(2),vsl( &
                     ldvsl,*),vsr(ldvsr,*),work(*)
           ! Function Arguments
           procedure(la_selctg_s) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl,wantsb, &
                     wantse,wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,ip,iright, &
                     irows,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -16
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -18
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'SGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'SORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'SORGQR',' ',n,1,n, &
                              -1))
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 6
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -22
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           safmin = la_slamch('S')
           safmax = one/safmin
           call la_slabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n for permutation parameters)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_sggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_slaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_sorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_slaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_sgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_shgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 50
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           ! (workspace: if ijob >= 1, need max( 8*(n+1), 2*sdim*(n-sdim) )
                       ! otherwise, need 8*(n+1) )
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              call la_stgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
              beta,vsl,ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork, &
                        liwork,ierr)
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -22) then
                  ! not enough real workspace
                 info = -22
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_sggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_sggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_slascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           50 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_sggesx
     !> DGGESX: computes for a pair of N-by-N real nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the real Schur form (S,T), and,
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_dggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim, &
     alphar,alphai,beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,iwork,liwork, &
               bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),rconde(2),rcondv(2),vsl( &
                     ldvsl,*),vsr(ldvsr,*),work(*)
           ! Function Arguments
           procedure(la_selctg_d) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl,wantsb, &
                     wantse,wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,ip,iright, &
                     irows,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -16
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -18
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'DGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'DORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'DORGQR',' ',n,1,n, &
                              -1))
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 6
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -22
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           safmin = la_dlamch('S')
           safmax = one/safmin
           call la_dlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n for permutation parameters)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_dggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_dlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_dorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_dlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_dgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_dhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 60
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           ! (workspace: if ijob >= 1, need max( 8*(n+1), 2*sdim*(n-sdim) )
                       ! otherwise, need 8*(n+1) )
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              call la_dtgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
              beta,vsl,ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork, &
                        liwork,ierr)
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -22) then
                  ! not enough real workspace
                 info = -22
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_dggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_dggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_dlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           60 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_dggesx
     !> QGGESX: computes for a pair of N-by-N real nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the real Schur form (S,T), and,
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_qggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim, &
     alphar,alphai,beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,iwork,liwork, &
               bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),rconde(2),rcondv(2),vsl( &
                     ldvsl,*),vsr(ldvsr,*),work(*)
           ! Function Arguments
           procedure(la_selctg_q) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl,wantsb, &
                     wantse,wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,ip,iright, &
                     irows,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -16
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -18
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = max(8*n,6*n + 16)
                 maxwrk = minwrk - n + n*la_ilaenv(1,'QGEQRF',' ',n,1,n,0)
                 maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'QORMQR',' ',n,1,n,-1 &
                           ))
                 if (ilvsl) then
                    maxwrk = max(maxwrk,minwrk - n + n*la_ilaenv(1,'QORGQR',' ',n,1,n, &
                              -1))
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 6
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -22
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           safmin = la_qlamch('S')
           safmax = one/safmin
           call la_qlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need 6*n + 2*n for permutation parameters)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_qggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_qlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_qorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_qlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_qgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (workspace: need n)
           iwrk = itau
           call la_qhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 60
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           ! (workspace: if ijob >= 1, need max( 8*(n+1), 2*sdim*(n-sdim) )
                       ! otherwise, need 8*(n+1) )
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              call la_qtgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
              beta,vsl,ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork, &
                        liwork,ierr)
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -22) then
                  ! not enough real workspace
                 info = -22
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_qggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_qggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_qlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           60 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_qggesx

     !> SGGEV: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_sggev(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
               ldvr,work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,maxwrk,minwrk
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              minwrk = max(1,8*n)
              maxwrk = max(1,n*(7 + la_ilaenv(1,'SGEQRF',' ',n,1,n,0)))
              maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'SORMQR',' ',n,1,n,0)))

              if (ilvl) then
                 maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'SORGQR',' ',n,1,n,-1)))

              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -16
           end if
           if (info /= 0) then
              call la_xerbla('SGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (workspace: need 6*n)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_sggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_slaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_sorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_slaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_sgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_sgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_shgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           ! (workspace: need 6*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_stgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_sggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_sggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_sggev
     !> DGGEV: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_dggev(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
               ldvr,work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,maxwrk,minwrk
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              minwrk = max(1,8*n)
              maxwrk = max(1,n*(7 + la_ilaenv(1,'DGEQRF',' ',n,1,n,0)))
              maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'DORMQR',' ',n,1,n,0)))

              if (ilvl) then
                 maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'DORGQR',' ',n,1,n,-1)))

              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (workspace: need 6*n)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_dggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_dlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_dorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_dlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_dgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_dgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_dhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           ! (workspace: need 6*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_dtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_dggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_dggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_dggev
     !> QGGEV: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_qggev(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
               ldvr,work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,maxwrk,minwrk
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              minwrk = max(1,8*n)
              maxwrk = max(1,n*(7 + la_ilaenv(1,'QGEQRF',' ',n,1,n,0)))
              maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'QORMQR',' ',n,1,n,0)))

              if (ilvl) then
                 maxwrk = max(maxwrk,n*(7 + la_ilaenv(1,'QORGQR',' ',n,1,n,-1)))

              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (workspace: need 6*n)
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_qggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (workspace: need n, prefer n*nb)
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_qlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_qorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_qlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_qgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_qgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_qhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           ! (workspace: need 6*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_qtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_qggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_qggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_qggev

     !> SGGEVX: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_sggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alphar,alphai, &
     beta,vl,ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork, &
               iwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(sp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),lscale(*),rconde(*),rcondv(*) &
                     ,rscale(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,pair,wantsb,wantse, &
                     wantsn,wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk,mm
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (noscl .or. la_lsame(balanc,'S') .or. la_lsame(balanc,'B'))) &
                     then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -14
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 if (noscl .and. .not. ilv) then
                    minwrk = 2*n
                 else
                    minwrk = 6*n
                 end if
                 if (wantse) then
                    minwrk = 10*n
                 else if (wantsv .or. wantsb) then
                    minwrk = 2*n*(n + 4) + 16
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'SGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'SORMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'SORGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -26
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_sggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,work,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_slange('1',n,n,a,lda,work(1))
           if (ilascl) then
              work(1) = abnrm
              call la_slascl('G',0,0,anrmto,anrm,1,1,work(1),1,ierr)
              abnrm = work(1)
           end if
           bbnrm = la_slange('1',n,n,b,ldb,work(1))
           if (ilbscl) then
              work(1) = bbnrm
              call la_slascl('G',0,0,bnrmto,bnrm,1,1,work(1),1,ierr)
              bbnrm = work(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to a
           ! (workspace: need n, prefer n*nb)
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_slaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_sorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_slaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_sgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_sgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_shgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work,lwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 130
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! (workspace: la_stgevc: need 6*n
                       ! la_stgsna: need 2*n*(n+2)+16 if sense = 'v' or 'b',
                               ! need n otherwise )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_stgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 130
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_stgevc) and estimate condition
                 ! numbers (la_stgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to recalculate
                 ! eigenvectors and estimate one condition numbers at a time.
                 pair = .false.
                 loop_20: do i = 1,n
                    if (pair) then
                       pair = .false.
                       cycle loop_20
                    end if
                    mm = 1
                    if (i < n) then
                       if (a(i + 1,i) /= zero) then
                          pair = .true.
                          mm = 2
                       end if
                    end if
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    if (mm == 1) then
                       bwork(i) = .true.
                    else if (mm == 2) then
                       bwork(i) = .true.
                       bwork(i + 1) = .true.
                    end if
                    iwrk = mm*n + 1
                    iwrk1 = iwrk + mm*n
                    ! compute a pair of left and right eigenvectors.
                    ! (compute workspace: need up to 4*n + 6*n)
                    if (wantse .or. wantsb) then
                       call la_stgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,mm,m,work(iwrk1),ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 130
                       end if
                    end if
                    call la_stgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),mm,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                               ierr)
                 end do loop_20
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_sggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_70: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_70
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_70
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                       vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_70
           end if
           if (ilvr) then
              call la_sggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_120: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_120
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_120
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                       vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_120
           end if
           ! undo scaling if necessary
           130 continue
           if (ilascl) then
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_sggevx
     !> DGGEVX: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_dggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alphar,alphai, &
     beta,vl,ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork, &
               iwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(dp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),lscale(*),rconde(*),rcondv(*) &
                     ,rscale(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,pair,wantsb,wantse, &
                     wantsn,wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk,mm
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') .or. &
                     la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -14
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 if (noscl .and. .not. ilv) then
                    minwrk = 2*n
                 else
                    minwrk = 6*n
                 end if
                 if (wantse .or. wantsb) then
                    minwrk = 10*n
                 end if
                 if (wantsv .or. wantsb) then
                    minwrk = max(minwrk,2*n*(n + 4) + 16)
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'DGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'DORMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'DORGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -26
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_dggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,work,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_dlange('1',n,n,a,lda,work(1))
           if (ilascl) then
              work(1) = abnrm
              call la_dlascl('G',0,0,anrmto,anrm,1,1,work(1),1,ierr)
              abnrm = work(1)
           end if
           bbnrm = la_dlange('1',n,n,b,ldb,work(1))
           if (ilbscl) then
              work(1) = bbnrm
              call la_dlascl('G',0,0,bnrmto,bnrm,1,1,work(1),1,ierr)
              bbnrm = work(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to a
           ! (workspace: need n, prefer n*nb)
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_dlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_dorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_dlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_dgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_dgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_dhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work,lwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 130
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! (workspace: la_dtgevc: need 6*n
                       ! la_dtgsna: need 2*n*(n+2)+16 if sense = 'v' or 'b',
                               ! need n otherwise )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_dtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 130
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_dtgevc) and estimate condition
                 ! numbers (la_dtgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to recalculate
                 ! eigenvectors and estimate one condition numbers at a time.
                 pair = .false.
                 loop_20: do i = 1,n
                    if (pair) then
                       pair = .false.
                       cycle loop_20
                    end if
                    mm = 1
                    if (i < n) then
                       if (a(i + 1,i) /= zero) then
                          pair = .true.
                          mm = 2
                       end if
                    end if
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    if (mm == 1) then
                       bwork(i) = .true.
                    else if (mm == 2) then
                       bwork(i) = .true.
                       bwork(i + 1) = .true.
                    end if
                    iwrk = mm*n + 1
                    iwrk1 = iwrk + mm*n
                    ! compute a pair of left and right eigenvectors.
                    ! (compute workspace: need up to 4*n + 6*n)
                    if (wantse .or. wantsb) then
                       call la_dtgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,mm,m,work(iwrk1),ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 130
                       end if
                    end if
                    call la_dtgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),mm,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                               ierr)
                 end do loop_20
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_dggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_70: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_70
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_70
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                       vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_70
           end if
           if (ilvr) then
              call la_dggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_120: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_120
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_120
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                       vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_120
           end if
           ! undo scaling if necessary
           130 continue
           if (ilascl) then
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_dggevx
     !> QGGEVX: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_qggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alphar,alphai, &
     beta,vl,ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork, &
               iwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(qp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),lscale(*),rconde(*),rcondv(*) &
                     ,rscale(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,pair,wantsb,wantse, &
                     wantsn,wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk,mm
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') .or. &
                     la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -14
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 if (noscl .and. .not. ilv) then
                    minwrk = 2*n
                 else
                    minwrk = 6*n
                 end if
                 if (wantse .or. wantsb) then
                    minwrk = 10*n
                 end if
                 if (wantsv .or. wantsb) then
                    minwrk = max(minwrk,2*n*(n + 4) + 16)
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'QGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'QORMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'QORGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -26
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_qggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,work,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_qlange('1',n,n,a,lda,work(1))
           if (ilascl) then
              work(1) = abnrm
              call la_qlascl('G',0,0,anrmto,anrm,1,1,work(1),1,ierr)
              abnrm = work(1)
           end if
           bbnrm = la_qlange('1',n,n,b,ldb,work(1))
           if (ilbscl) then
              work(1) = bbnrm
              call la_qlascl('G',0,0,bnrmto,bnrm,1,1,work(1),1,ierr)
              bbnrm = work(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to a
           ! (workspace: need n, prefer n*nb)
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_qlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_qorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_qlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_qgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_qgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (workspace: need n)
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_qhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work,lwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 130
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! (workspace: la_qtgevc: need 6*n
                       ! la_qtgsna: need 2*n*(n+2)+16 if sense = 'v' or 'b',
                               ! need n otherwise )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_qtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 130
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_qtgevc) and estimate condition
                 ! numbers (la_qtgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to recalculate
                 ! eigenvectors and estimate one condition numbers at a time.
                 pair = .false.
                 loop_20: do i = 1,n
                    if (pair) then
                       pair = .false.
                       cycle loop_20
                    end if
                    mm = 1
                    if (i < n) then
                       if (a(i + 1,i) /= zero) then
                          pair = .true.
                          mm = 2
                       end if
                    end if
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    if (mm == 1) then
                       bwork(i) = .true.
                    else if (mm == 2) then
                       bwork(i) = .true.
                       bwork(i + 1) = .true.
                    end if
                    iwrk = mm*n + 1
                    iwrk1 = iwrk + mm*n
                    ! compute a pair of left and right eigenvectors.
                    ! (compute workspace: need up to 4*n + 6*n)
                    if (wantse .or. wantsb) then
                       call la_qtgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,mm,m,work(iwrk1),ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 130
                       end if
                    end if
                    call la_qtgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),mm,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                               ierr)
                 end do loop_20
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_qggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_70: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_70
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_70
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                       vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_70
           end if
           if (ilvr) then
              call la_qggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_120: do jc = 1,n
                 if (alphai(jc) < zero) cycle loop_120
                 temp = zero
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)))
                    end do
                 else
                    do jr = 1,n
                       temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                    end do
                 end if
                 if (temp < smlnum) cycle loop_120
                 temp = one/temp
                 if (alphai(jc) == zero) then
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 else
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                       vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                    end do
                 end if
              end do loop_120
           end if
           ! undo scaling if necessary
           130 continue
           if (ilascl) then
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = maxwrk
           return
     end subroutine la_qggevx

     !> SGEES: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A matrix is in real Schur form if it is upper quasi-triangular with
     !> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
     !> form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_sgees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork, &
               bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_s) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,maxwrk,minwrk
           real(sp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_shseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'SGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_shseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'SORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_slascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need n)
           ibal = 1
           call la_sgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_sgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_slacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (workspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_shseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_slascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_slascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (workspace: none needed)
              call la_strsen('N',jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,s,sep, &
                        work(iwrk),lwork - iwrk + 1,idum,1,icond)
              if (icond > 0) info = n + icond
           end if
           if (wantvs) then
              ! undo balancing
              ! (workspace: need n)
              call la_sgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_slascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_scopy(n,a,lda + 1,wr,1)
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wi,max(ilo - 1,1), &
                              ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_sswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_sswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                             call la_sswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              ! undo scaling for the imaginary part of the eigenvalues
              call la_slascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           return
     end subroutine la_sgees
     !> DGEES: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A matrix is in real Schur form if it is upper quasi-triangular with
     !> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
     !> form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_dgees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork, &
               bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_d) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,maxwrk,minwrk
           real(dp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_dhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'DGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_dhseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'DORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_dlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need n)
           ibal = 1
           call la_dgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_dgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_dlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (workspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_dhseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_dlascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_dlascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (workspace: none needed)
              call la_dtrsen('N',jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,s,sep, &
                        work(iwrk),lwork - iwrk + 1,idum,1,icond)
              if (icond > 0) info = n + icond
           end if
           if (wantvs) then
              ! undo balancing
              ! (workspace: need n)
              call la_dgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_dlascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_dcopy(n,a,lda + 1,wr,1)
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,max(ilo - 1,1), &
                              ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_dswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_dswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                             call la_dswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              ! undo scaling for the imaginary part of the eigenvalues
              call la_dlascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           return
     end subroutine la_dgees
     !> QGEES: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A matrix is in real Schur form if it is upper quasi-triangular with
     !> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
     !> form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_qgees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,lwork, &
               bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_q) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,maxwrk,minwrk
           real(qp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_qhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'QGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_qhseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'QORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_qlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (workspace: need n)
           ibal = 1
           call la_qgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_qgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_qlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (workspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_qhseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_qlascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_qlascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (workspace: none needed)
              call la_qtrsen('N',jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,s,sep, &
                        work(iwrk),lwork - iwrk + 1,idum,1,icond)
              if (icond > 0) info = n + icond
           end if
           if (wantvs) then
              ! undo balancing
              ! (workspace: need n)
              call la_qgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_qlascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_qcopy(n,a,lda + 1,wr,1)
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,max(ilo - 1,1), &
                              ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_qswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_qswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                             call la_qswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              ! undo scaling for the imaginary part of the eigenvalues
              call la_qlascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           return
     end subroutine la_qgees

     !> SGEESX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_sp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A real matrix is in real Schur form if it is upper quasi-triangular
     !> with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
     !> the form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_sgeesx(jobvs,sort,select,sense,n,a,lda,sdim,wr,wi,vs,ldvs, &
               rconde,rcondv,work,lwork,iwork,liwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,liwork,lwork,n
           real(sp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_s) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantsb,wantse,wantsn,wantst, &
                     wantsv,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,lwrk,liwrk,maxwrk,minwrk
           real(sp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "rworkspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! iworkspace refers to integer workspace.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_shseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_strsen later
             ! in the code.)
           if (info == 0) then
              liwrk = 1
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'SGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_shseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'SORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk,n + (n*n)/2)
                 if (wantsv .or. wantsb) liwrk = (n*n)/4
              end if
              iwork(1) = liwrk
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -16
              else if (liwork < 1 .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_slascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (rworkspace: need n)
           ibal = 1
           call la_sgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (rworkspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_sgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_slacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (rworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (rworkspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_shseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_slascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_slascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (rworkspace: if sense is not 'n', need n+2*sdim*(n-sdim)
                           ! otherwise, need n )
              ! (iworkspace: if sense is 'v' or 'b', need sdim*(n-sdim)
                           ! otherwise, need 0 )
              call la_strsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,iwork,liwork,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,n + 2*sdim*(n - sdim))
              if (icond == -15) then
                 ! not enough real workspace
                 info = -16
              else if (icond == -17) then
                 ! not enough integer workspace
                 info = -18
              else if (icond > 0) then
                 ! la_strsen failed to reorder or to restore standard schur form
                 info = icond + n
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (rworkspace: need n)
              call la_sgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_slascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_scopy(n,a,lda + 1,wr,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_slascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_sswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_sswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                            call la_sswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              call la_slascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           if (wantsv .or. wantsb) then
              iwork(1) = sdim*(n - sdim)
           else
              iwork(1) = 1
           end if
           return
     end subroutine la_sgeesx
     !> DGEESX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_dp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A real matrix is in real Schur form if it is upper quasi-triangular
     !> with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
     !> the form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_dgeesx(jobvs,sort,select,sense,n,a,lda,sdim,wr,wi,vs,ldvs, &
               rconde,rcondv,work,lwork,iwork,liwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,liwork,lwork,n
           real(dp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_d) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantsb,wantse,wantsn,wantst, &
                     wantsv,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,liwrk,lwrk,maxwrk,minwrk
           real(dp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "rworkspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! iworkspace refers to integer workspace.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_dhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_dtrsen later
             ! in the code.)
           if (info == 0) then
              liwrk = 1
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'DGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_dhseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'DORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk,n + (n*n)/2)
                 if (wantsv .or. wantsb) liwrk = (n*n)/4
              end if
              iwork(1) = liwrk
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -16
              else if (liwork < 1 .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_dlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (rworkspace: need n)
           ibal = 1
           call la_dgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (rworkspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_dgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_dlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (rworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (rworkspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_dhseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_dlascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_dlascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (rworkspace: if sense is not 'n', need n+2*sdim*(n-sdim)
                           ! otherwise, need n )
              ! (iworkspace: if sense is 'v' or 'b', need sdim*(n-sdim)
                           ! otherwise, need 0 )
              call la_dtrsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,iwork,liwork,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,n + 2*sdim*(n - sdim))
              if (icond == -15) then
                 ! not enough real workspace
                 info = -16
              else if (icond == -17) then
                 ! not enough integer workspace
                 info = -18
              else if (icond > 0) then
                 ! la_dtrsen failed to reorder or to restore standard schur form
                 info = icond + n
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (rworkspace: need n)
              call la_dgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_dlascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_dcopy(n,a,lda + 1,wr,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_dlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_dswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_dswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                            call la_dswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              call la_dlascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           if (wantsv .or. wantsb) then
              iwork(1) = max(1,sdim*(n - sdim))
           else
              iwork(1) = 1
           end if
           return
     end subroutine la_dgeesx
     !> QGEESX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues, the real Schur form T, and, optionally, the matrix of
     !> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> real Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_qp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A real matrix is in real Schur form if it is upper quasi-triangular
     !> with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
     !> the form
     !> [  a  b  ]
     !> [  c  a  ]
     !> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).

     subroutine la_qgeesx(jobvs,sort,select,sense,n,a,lda,sdim,wr,wi,vs,ldvs, &
               rconde,rcondv,work,lwork,iwork,liwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,liwork,lwork,n
           real(qp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: vs(ldvs,*),wi(*),work(*),wr(*)
           ! Function Arguments
           procedure(la_select_q) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,lastsl,lquery,lst2sl,scalea,wantsb,wantse,wantsn,wantst, &
                     wantsv,wantvs
           integer(ilp) :: hswork,i,i1,i2,ibal,icond,ierr,ieval,ihi,ilo,inxt,ip,itau, &
                     iwrk,liwrk,lwrk,maxwrk,minwrk
           real(qp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "rworkspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! iworkspace refers to integer workspace.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_qhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_qtrsen later
             ! in the code.)
           if (info == 0) then
              liwrk = 1
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'QGEHRD',' ',n,1,n,0)
                 minwrk = 3*n
                 call la_qhseqr('S',jobvs,n,1,n,a,lda,wr,wi,vs,ldvs,work,-1, &
                           ieval)
                 hswork = work(1)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,n + hswork)
                 else
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'QORGHR',' ',n,1,n, &
                               -1))
                    maxwrk = max(maxwrk,n + hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk,n + (n*n)/2)
                 if (wantsv .or. wantsb) liwrk = (n*n)/4
              end if
              iwork(1) = liwrk
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -16
              else if (liwork < 1 .and. .not. lquery) then
                 info = -18
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_qlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (rworkspace: need n)
           ibal = 1
           call la_qgebal('P',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (rworkspace: need 3*n, prefer 2*n+n*nb)
           itau = n + ibal
           iwrk = n + itau
           call la_qgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_qlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate orthogonal matrix in vs
              ! (rworkspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (rworkspace: need n+1, prefer n+hswork (see comments) )
           iwrk = itau
           call la_qhseqr('S',jobvs,n,ilo,ihi,a,lda,wr,wi,vs,ldvs,work(iwrk), &
                     lwork - iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) then
                 call la_qlascl('G',0,0,cscale,anrm,n,1,wr,n,ierr)
                 call la_qlascl('G',0,0,cscale,anrm,n,1,wi,n,ierr)
              end if
              do i = 1,n
                 bwork(i) = select(wr(i),wi(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (rworkspace: if sense is not 'n', need n+2*sdim*(n-sdim)
                           ! otherwise, need n )
              ! (iworkspace: if sense is 'v' or 'b', need sdim*(n-sdim)
                           ! otherwise, need 0 )
              call la_qtrsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,wr,wi,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,iwork,liwork,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,n + 2*sdim*(n - sdim))
              if (icond == -15) then
                 ! not enough real workspace
                 info = -16
              else if (icond == -17) then
                 ! not enough integer workspace
                 info = -18
              else if (icond > 0) then
                 ! la_qtrsen failed to reorder or to restore standard schur form
                 info = icond + n
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (rworkspace: need n)
              call la_qgebak('P','R',n,ilo,ihi,work(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_qlascl('H',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_qcopy(n,a,lda + 1,wr,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_qlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
              if (cscale == smlnum) then
                 ! if scaling back towards underflow, adjust wi if an
                 ! offdiagonal element of a 2-by-2 block in the schur form
                 ! underflows.
                 if (ieval > 0) then
                    i1 = ieval + 1
                    i2 = ihi - 1
                    call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
                 else if (wantst) then
                    i1 = 1
                    i2 = n - 1
                 else
                    i1 = ilo
                    i2 = ihi - 1
                 end if
                 inxt = i1 - 1
                 loop_20: do i = i1,i2
                    if (i < inxt) cycle loop_20
                    if (wi(i) == zero) then
                       inxt = i + 1
                    else
                       if (a(i + 1,i) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                       else if (a(i + 1,i) /= zero .and. a(i,i + 1) == zero) then
                          wi(i) = zero
                          wi(i + 1) = zero
                          if (i > 1) call la_qswap(i - 1,a(1,i),1,a(1,i + 1),1)
                          if (n > i + 1) call la_qswap(n - i - 1,a(i,i + 2),lda,a(i + 1,i + 2), &
                                    lda)
                          if (wantvs) then
                            call la_qswap(n,vs(1,i),1,vs(1,i + 1),1)
                          end if
                          a(i,i + 1) = a(i + 1,i)
                          a(i + 1,i) = zero
                       end if
                       inxt = i + 2
                    end if
                 end do loop_20
              end if
              call la_qlascl('G',0,0,cscale,anrm,n - ieval,1,wi(ieval + 1),max(n - ieval, &
                         1),ierr)
           end if
           if (wantst .and. info == 0) then
              ! check if reordering successful
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = select(wr(i),wi(i))
                 if (wi(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           work(1) = maxwrk
           if (wantsv .or. wantsb) then
              iwork(1) = max(1,sdim*(n - sdim))
           else
              iwork(1) = 1
           end if
           return
     end subroutine la_qgeesx

     !> SGEEV: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_sgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork, &
               info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: vl(ldvl,*),vr(ldvr,*),wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,itau,iwrk,k,lwork_trevc,maxwrk, &
                     minwrk,nout
           real(sp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -9
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_shseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'SGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'SORGHR',' ',n,1,n, &
                               -1))
                    call la_shseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_strevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else if (wantvr) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'SORGHR',' ',n,1,n, &
                               -1))
                    call la_shseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_strevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else
                    minwrk = 3*n
                    call la_shseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_slascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (workspace: need n)
           ibal = 1
           call la_sgebal('B',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = ibal + n
           iwrk = itau + n
           call la_sgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_slacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_shseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_slacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_slacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_shseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_shseqr('E','N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_shseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 4*n, prefer n + n + 2*n*nb)
              call la_strevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (workspace: need n)
              call la_sgebak('B','L',n,ilo,ihi,work(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_snrm2(n,vl(1,i),1)
                    call la_sscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_slapy2(la_snrm2(n,vl(1,i),1),la_snrm2(n, &
                              vl(1,i + 1),1))
                    call la_sscal(n,scl,vl(1,i),1)
                    call la_sscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_isamax(n,work(iwrk),1)
                    call la_slartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_srot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (workspace: need n)
              call la_sgebak('B','R',n,ilo,ihi,work(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_snrm2(n,vr(1,i),1)
                    call la_sscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_slapy2(la_snrm2(n,vr(1,i),1),la_snrm2(n, &
                              vr(1,i + 1),1))
                    call la_sscal(n,scl,vr(1,i),1)
                    call la_sscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_isamax(n,work(iwrk),1)
                    call la_slartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_srot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_slascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_slascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info > 0) then
                 call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_sgeev
     !> DGEEV: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork, &
               info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: vl(ldvl,*),vr(ldvr,*),wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,itau,iwrk,k,lwork_trevc,maxwrk, &
                     minwrk,nout
           real(dp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -9
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_dhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'DGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'DORGHR',' ',n,1,n, &
                               -1))
                    call la_dhseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_dtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else if (wantvr) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'DORGHR',' ',n,1,n, &
                               -1))
                    call la_dhseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_dtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else
                    minwrk = 3*n
                    call la_dhseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_dlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (workspace: need n)
           ibal = 1
           call la_dgebal('B',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = ibal + n
           iwrk = itau + n
           call la_dgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_dlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_dhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_dlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_dlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_dhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_dhseqr('E','N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_dhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 4*n, prefer n + n + 2*n*nb)
              call la_dtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (workspace: need n)
              call la_dgebak('B','L',n,ilo,ihi,work(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_dnrm2(n,vl(1,i),1)
                    call la_dscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_dlapy2(la_dnrm2(n,vl(1,i),1),la_dnrm2(n, &
                              vl(1,i + 1),1))
                    call la_dscal(n,scl,vl(1,i),1)
                    call la_dscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_idamax(n,work(iwrk),1)
                    call la_dlartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_drot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (workspace: need n)
              call la_dgebak('B','R',n,ilo,ihi,work(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_dnrm2(n,vr(1,i),1)
                    call la_dscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_dlapy2(la_dnrm2(n,vr(1,i),1),la_dnrm2(n, &
                              vr(1,i + 1),1))
                    call la_dscal(n,scl,vr(1,i),1)
                    call la_dscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_idamax(n,work(iwrk),1)
                    call la_dlartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_drot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_dlascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_dlascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info > 0) then
                 call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_dgeev
     !> QGEEV: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_qgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork, &
               info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: vl(ldvl,*),vr(ldvr,*),wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,itau,iwrk,k,lwork_trevc,maxwrk, &
                     minwrk,nout
           real(qp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -9
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_qhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = 2*n + n*la_ilaenv(1,'QGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'QORGHR',' ',n,1,n, &
                               -1))
                    call la_qhseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_qtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else if (wantvr) then
                    minwrk = 4*n
                    maxwrk = max(maxwrk,2*n + (n - 1)*la_ilaenv(1,'QORGHR',' ',n,1,n, &
                               -1))
                    call la_qhseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                    call la_qtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    maxwrk = max(maxwrk,4*n)
                 else
                    minwrk = 3*n
                    call la_qhseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                    hswork = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + 1,n + hswork)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -13
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_qlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (workspace: need n)
           ibal = 1
           call la_qgebal('B',n,a,lda,ilo,ihi,work(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (workspace: need 3*n, prefer 2*n+n*nb)
           itau = ibal + n
           iwrk = itau + n
           call la_qgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_qlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_qhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_qlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_qlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 3*n-1, prefer 2*n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_qhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (workspace: need n+1, prefer n+hswork (see comments) )
              iwrk = itau
              call la_qhseqr('E','N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_qhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 4*n, prefer n + n + 2*n*nb)
              call la_qtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (workspace: need n)
              call la_qgebak('B','L',n,ilo,ihi,work(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_qnrm2(n,vl(1,i),1)
                    call la_qscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_qlapy2(la_qnrm2(n,vl(1,i),1),la_qnrm2(n, &
                              vl(1,i + 1),1))
                    call la_qscal(n,scl,vl(1,i),1)
                    call la_qscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_iqamax(n,work(iwrk),1)
                    call la_qlartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_qrot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (workspace: need n)
              call la_qgebak('B','R',n,ilo,ihi,work(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_qnrm2(n,vr(1,i),1)
                    call la_qscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_qlapy2(la_qnrm2(n,vr(1,i),1),la_qnrm2(n, &
                              vr(1,i + 1),1))
                    call la_qscal(n,scl,vr(1,i),1)
                    call la_qscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(iwrk + k - 1) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_iqamax(n,work(iwrk),1)
                    call la_qlartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_qrot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_qlascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_qlascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info > 0) then
                 call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_qgeev

     !> SGEEVX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_sp of the LAPACK
     !> Users' Guide.

     subroutine la_sgeevx(balanc,jobvl,jobvr,sense,n,a,lda,wr,wi,vl,ldvl,vr,ldvr, &
               ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(sp),intent(out) :: abnrm
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: rconde(*),rcondv(*),scale(*),vl(ldvl,*),vr(ldvr,*),wi(*), &
                      work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(sp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') .or. &
                     la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_shseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'SGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_strevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_shseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                 else if (wantvr) then
                    call la_strevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_shseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                 else
                    if (wntsnn) then
                       call la_shseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    else
                       call la_shseqr('S','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. wntsnn) minwrk = max(minwrk,n*n + 6*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. wntsnn) maxwrk = max(maxwrk,n*n + 6*n)
                 else
                    minwrk = 3*n
                    if ((.not. wntsnn) .and. (.not. wntsne)) minwrk = max(minwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'SORGHR',' ',n,1,n,- &
                              1))
                    if ((.not. wntsnn) .and. (.not. wntsne)) maxwrk = max(maxwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,3*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_slange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_slascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_sgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_slange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_slascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (workspace: need 2*n, prefer n+n*nb)
           itau = 1
           iwrk = itau + n
           call la_sgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_slacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_shseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_slacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_slacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_sorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_shseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_shseqr(job,'N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_shseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 3*n, prefer n + 2*n*nb)
              call la_strevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           ! compute condition numbers if desired
           ! (workspace: need n*n+6*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_strsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,iwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_sgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_snrm2(n,vl(1,i),1)
                    call la_sscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_slapy2(la_snrm2(n,vl(1,i),1),la_snrm2(n, &
                              vl(1,i + 1),1))
                    call la_sscal(n,scl,vl(1,i),1)
                    call la_sscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(k) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_isamax(n,work,1)
                    call la_slartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_srot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_sgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_snrm2(n,vr(1,i),1)
                    call la_sscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_slapy2(la_snrm2(n,vr(1,i),1),la_snrm2(n, &
                              vr(1,i + 1),1))
                    call la_sscal(n,scl,vr(1,i),1)
                    call la_sscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(k) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_isamax(n,work,1)
                    call la_slartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_srot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_slascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_slascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_slascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_slascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_sgeevx
     !> DGEEVX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_dp of the LAPACK
     !> Users' Guide.

     subroutine la_dgeevx(balanc,jobvl,jobvr,sense,n,a,lda,wr,wi,vl,ldvl,vr,ldvr, &
               ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(dp),intent(out) :: abnrm
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: rconde(*),rcondv(*),scale(*),vl(ldvl,*),vr(ldvr,*),wi(*), &
                      work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(dp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') .or. &
                     la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_dhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'DGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_dtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_dhseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                 else if (wantvr) then
                    call la_dtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_dhseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                 else
                    if (wntsnn) then
                       call la_dhseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    else
                       call la_dhseqr('S','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. wntsnn) minwrk = max(minwrk,n*n + 6*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. wntsnn) maxwrk = max(maxwrk,n*n + 6*n)
                 else
                    minwrk = 3*n
                    if ((.not. wntsnn) .and. (.not. wntsne)) minwrk = max(minwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'DORGHR',' ',n,1,n,- &
                              1))
                    if ((.not. wntsnn) .and. (.not. wntsne)) maxwrk = max(maxwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,3*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_dlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_dlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_dgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_dlange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_dlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (workspace: need 2*n, prefer n+n*nb)
           itau = 1
           iwrk = itau + n
           call la_dgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_dlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_dhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_dlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_dlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_dorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_dhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_dhseqr(job,'N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_dhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 3*n, prefer n + 2*n*nb)
              call la_dtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           ! compute condition numbers if desired
           ! (workspace: need n*n+6*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_dtrsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,iwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_dgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_dnrm2(n,vl(1,i),1)
                    call la_dscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_dlapy2(la_dnrm2(n,vl(1,i),1),la_dnrm2(n, &
                              vl(1,i + 1),1))
                    call la_dscal(n,scl,vl(1,i),1)
                    call la_dscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(k) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_idamax(n,work,1)
                    call la_dlartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_drot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_dgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_dnrm2(n,vr(1,i),1)
                    call la_dscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_dlapy2(la_dnrm2(n,vr(1,i),1),la_dnrm2(n, &
                              vr(1,i + 1),1))
                    call la_dscal(n,scl,vr(1,i),1)
                    call la_dscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(k) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_idamax(n,work,1)
                    call la_dlartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_drot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_dlascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_dlascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_dlascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_dlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_dgeevx
     !> QGEEVX: computes for an N-by-N real nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate-transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_qp of the LAPACK
     !> Users' Guide.

     subroutine la_qgeevx(balanc,jobvl,jobvr,sense,n,a,lda,wr,wi,vl,ldvl,vr,ldvr, &
               ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(qp),intent(out) :: abnrm
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: rconde(*),rcondv(*),scale(*),vl(ldvl,*),vr(ldvr,*),wi(*), &
                      work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(qp) :: anrm,bignum,cs,cscale,eps,r,scl,smlnum,sn
           ! Local Arrays
           logical(lk) :: select(1)
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') .or. &
                     la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_qhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'QGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_qtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_qhseqr('S','V',n,1,n,a,lda,wr,wi,vl,ldvl,work,-1, &
                              info)
                 else if (wantvr) then
                    call la_qtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_qhseqr('S','V',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                              info)
                 else
                    if (wntsnn) then
                       call la_qhseqr('E','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    else
                       call la_qhseqr('S','N',n,1,n,a,lda,wr,wi,vr,ldvr,work,-1, &
                                 info)
                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. wntsnn) minwrk = max(minwrk,n*n + 6*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. wntsnn) maxwrk = max(maxwrk,n*n + 6*n)
                 else
                    minwrk = 3*n
                    if ((.not. wntsnn) .and. (.not. wntsne)) minwrk = max(minwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'QORGHR',' ',n,1,n,- &
                              1))
                    if ((.not. wntsnn) .and. (.not. wntsne)) maxwrk = max(maxwrk,n*n + 6*n)

                    maxwrk = max(maxwrk,3*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_qlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_qlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_qgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_qlange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_qlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (workspace: need 2*n, prefer n+n*nb)
           itau = 1
           iwrk = itau + n
           call la_qgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_qlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate orthogonal matrix in vl
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_qhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vl,ldvl,work(iwrk), &
                        lwork - iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_qlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_qlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate orthogonal matrix in vr
              ! (workspace: need 2*n-1, prefer n+(n-1)*nb)
              call la_qorghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_qhseqr('S','V',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (workspace: need 1, prefer hswork (see comments) )
              iwrk = itau
              call la_qhseqr(job,'N',n,ilo,ihi,a,lda,wr,wi,vr,ldvr,work(iwrk), &
                        lwork - iwrk + 1,info)
           end if
           ! if info /= 0 from la_qhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (workspace: need 3*n, prefer n + 2*n*nb)
              call la_qtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,ierr)
           end if
           ! compute condition numbers if desired
           ! (workspace: need n*n+6*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_qtrsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,iwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_qgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_qnrm2(n,vl(1,i),1)
                    call la_qscal(n,scl,vl(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_qlapy2(la_qnrm2(n,vl(1,i),1),la_qnrm2(n, &
                              vl(1,i + 1),1))
                    call la_qscal(n,scl,vl(1,i),1)
                    call la_qscal(n,scl,vl(1,i + 1),1)
                    do k = 1,n
                       work(k) = vl(k,i)**2 + vl(k,i + 1)**2
                    end do
                    k = la_iqamax(n,work,1)
                    call la_qlartg(vl(k,i),vl(k,i + 1),cs,sn,r)
                    call la_qrot(n,vl(1,i),1,vl(1,i + 1),1,cs,sn)
                    vl(k,i + 1) = zero
                 end if
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_qgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 if (wi(i) == zero) then
                    scl = one/la_qnrm2(n,vr(1,i),1)
                    call la_qscal(n,scl,vr(1,i),1)
                 else if (wi(i) > zero) then
                    scl = one/la_qlapy2(la_qnrm2(n,vr(1,i),1),la_qnrm2(n, &
                              vr(1,i + 1),1))
                    call la_qscal(n,scl,vr(1,i),1)
                    call la_qscal(n,scl,vr(1,i + 1),1)
                    do k = 1,n
                       work(k) = vr(k,i)**2 + vr(k,i + 1)**2
                    end do
                    k = la_iqamax(n,work,1)
                    call la_qlartg(vr(k,i),vr(k,i + 1),cs,sn,r)
                    call la_qrot(n,vr(1,i),1,vr(1,i + 1),1,cs,sn)
                    vr(k,i + 1) = zero
                 end if
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_qlascl('G',0,0,cscale,anrm,n - info,1,wr(info + 1),max(n - info,1 &
                        ),ierr)
              call la_qlascl('G',0,0,cscale,anrm,n - info,1,wi(info + 1),max(n - info,1 &
                        ),ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_qlascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wr,n,ierr)
                 call la_qlascl('G',0,0,cscale,anrm,ilo - 1,1,wi,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_qgeevx

     !> SGGES3: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> SGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_sgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_s) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           else if (lwork < 6*n + 16 .and. .not. lquery) then
              info = -19
           end if
           ! compute workspace
           if (info == 0) then
              call la_sgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(6*n + 16,3*n + int(work(1),KIND=ilp))
              call la_sormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_sorgqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              end if
              call la_sgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              call la_slaqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                        beta,vsl,ldvsl,vsr,ldvsr,work,-1,0,ierr)
              lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              if (wantst) then
                 call la_stgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
                 beta,vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)

                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           safmin = la_slamch('S')
           safmax = one/safmin
           call la_slabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_sggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_slaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_sorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_slaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_sgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_slaqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 40
           end if
           ! sort eigenvalues alpha/beta if desired
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_stgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_sggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_sggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_slascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           40 continue
           work(1) = lwkopt
           return
     end subroutine la_sgges3
     !> DGGES3: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> DGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_dgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_d) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           else if (lwork < 6*n + 16 .and. .not. lquery) then
              info = -19
           end if
           ! compute workspace
           if (info == 0) then
              call la_dgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(6*n + 16,3*n + int(work(1),KIND=ilp))
              call la_dormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_dorgqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              end if
              call la_dgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              call la_dlaqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                        beta,vsl,ldvsl,vsr,ldvsr,work,-1,0,ierr)
              lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              if (wantst) then
                 call la_dtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
                 beta,vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)

                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           safmin = la_dlamch('S')
           safmax = one/safmin
           call la_dlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_dggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_dlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_dorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_dlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_dgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_dlaqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 50
           end if
           ! sort eigenvalues alpha/beta if desired
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_dtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_dggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_dggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_dlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           50 continue
           work(1) = lwkopt
           return
     end subroutine la_dgges3
     !> QGGES3: computes for a pair of N-by-N real nonsymmetric matrices (A,B),
     !> the generalized eigenvalues, the generalized real Schur form (S,T),
     !> optionally, the left and/or right matrices of Schur vectors (VSL and
     !> VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> quasi-triangular matrix S and the upper triangular matrix T.The
     !> leading columns of VSL and VSR then form an orthonormal basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> QGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or both being zero.
     !> A pair of matrices (S,T) is in generalized real Schur form if T is
     !> upper triangular with non-negative diagonal and S is block upper
     !> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
     !> to real generalized eigenvalues, while 2-by-2 blocks of S will be
     !> "standardized" by making the corresponding elements of T have the
     !> form:
     !> [  a  0  ]
     !> [  0  b  ]
     !> and the pair of corresponding 2-by-2 blocks in S and T will have a
     !> complex conjugate pair of generalized eigenvalues.

     subroutine la_qgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alphar, &
               alphai,beta,vsl,ldvsl,vsr,ldvsr,work,lwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*), &
                     work(*)
           ! Function Arguments
           procedure(la_selctg_q) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,lst2sl, &
                     wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,ip,iright,irows, &
                     itau,iwrk,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,safmax,safmin, &
                     smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           else if (lwork < 6*n + 16 .and. .not. lquery) then
              info = -19
           end if
           ! compute workspace
           if (info == 0) then
              call la_qgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(6*n + 16,3*n + int(work(1),KIND=ilp))
              call la_qormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_qorgqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              end if
              call la_qgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              call la_qlaqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                        beta,vsl,ldvsl,vsr,ldvsr,work,-1,0,ierr)
              lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              if (wantst) then
                 call la_qtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai, &
                 beta,vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)

                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           safmin = la_qlamch('S')
           safmax = one/safmin
           call la_qlabad(safmin,safmax)
           smlnum = sqrt(safmin)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_qggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = iwrk
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_qlaset('FULL',n,n,zero,one,vsl,ldvsl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_qorgqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_qlaset('FULL',n,n,zero,one,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_qgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_qlaqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vsl,ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 50
           end if
           ! sort eigenvalues alpha/beta if desired
           sdim = 0
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) then
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
                 call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
              end if
              if (ilbscl) call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alphar(i),alphai(i),beta(i))
              end do
              call la_qtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alphar,alphai,beta, &
              vsl,ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1, &
                        ierr)
              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_qggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsl,ldvsl,ierr)
           if (ilvsr) call la_qggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n, &
                     vsr,ldvsr,ierr)
           ! check if unscaling would cause over/underflow, if so, rescale
           ! (alphar(i),alphai(i),beta(i)) so beta(i) is on the order of
           ! b(i,i) and alphar(i) and alphai(i) are on the order of a(i,i)
           if (ilascl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((alphar(i)/safmax) > (anrmto/anrm) .or. (safmin/alphar(i)) > ( &
                              anrm/anrmto)) then
                       work(1) = abs(a(i,i)/alphar(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    else if ((alphai(i)/safmax) > (anrmto/anrm) .or. (safmin/alphai(i) &
                               ) > (anrm/anrmto)) then
                       work(1) = abs(a(i,i + 1)/alphai(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           if (ilbscl) then
              do i = 1,n
                 if (alphai(i) /= zero) then
                    if ((beta(i)/safmax) > (bnrmto/bnrm) .or. (safmin/beta(i)) > ( &
                              bnrm/bnrmto)) then
                       work(1) = abs(b(i,i)/beta(i))
                       beta(i) = beta(i)*work(1)
                       alphar(i) = alphar(i)*work(1)
                       alphai(i) = alphai(i)*work(1)
                    end if
                 end if
              end do
           end if
           ! undo scaling
           if (ilascl) then
              call la_qlascl('H',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              lst2sl = .true.
              sdim = 0
              ip = 0
              do i = 1,n
                 cursl = selctg(alphar(i),alphai(i),beta(i))
                 if (alphai(i) == zero) then
                    if (cursl) sdim = sdim + 1
                    ip = 0
                    if (cursl .and. .not. lastsl) info = n + 2
                 else
                    if (ip == 1) then
                       ! last eigenvalue of conjugate pair
                       cursl = cursl .or. lastsl
                       lastsl = cursl
                       if (cursl) sdim = sdim + 2
                       ip = -1
                       if (cursl .and. .not. lst2sl) info = n + 2
                    else
                       ! first eigenvalue of conjugate pair
                       ip = 1
                    end if
                 end if
                 lst2sl = lastsl
                 lastsl = cursl
              end do
           end if
           50 continue
           work(1) = lwkopt
           return
     end subroutine la_qgges3

     !> SGGEV3: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_sggev3(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
                ldvr,work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           else if (lwork < max(1,8*n) .and. .not. lquery) then
              info = -16
           end if
           ! compute workspace
           if (info == 0) then
              call la_sgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,8*n,3*n + int(work(1),KIND=ilp))
              call la_sormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              call la_sgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,work, &
                        -1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_sorgqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
                 call la_slaqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              else
                 call la_slaqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = real(lwkopt,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('SGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_slange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_slascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_slange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_slascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_sggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_sgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_sormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_slaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_slacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_sorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_slaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_sgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_sgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_slaqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_stgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_sggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_sggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_slascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_slascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = real(lwkopt,KIND=sp)
           return
     end subroutine la_sggev3
     !> DGGEV3: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_dggev3(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
                ldvr,work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           else if (lwork < max(1,8*n) .and. .not. lquery) then
              info = -16
           end if
           ! compute workspace
           if (info == 0) then
              call la_dgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,8*n,3*n + int(work(1),KIND=ilp))
              call la_dormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_dorgqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              end if
              if (ilv) then
                 call la_dgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
                 call la_dlaqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              else
                 call la_dgghd3('N','N',n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,work,- &
                           1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
                 call la_dlaqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_dlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_dlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_dlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_dlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_dggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_dgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_dormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_dlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_dlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_dorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_dlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_dgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_dgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_dlaqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_dtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_dggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_dggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_dlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_dlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = lwkopt
           return
     end subroutine la_dggev3
     !> QGGEV3: computes for a pair of N-by-N real nonsymmetric matrices (A,B)
     !> the generalized eigenvalues, and optionally, the left and/or right
     !> generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B .
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_qggev3(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,vr, &
                ldvr,work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: alphai(*),alphar(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)

        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,itau, &
                     iwrk,jc,jr,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -12
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -14
           else if (lwork < max(1,8*n) .and. .not. lquery) then
              info = -16
           end if
           ! compute workspace
           if (info == 0) then
              call la_qgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,8*n,3*n + int(work(1),KIND=ilp))
              call la_qormqr('L','T',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_qorgqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
              end if
              if (ilv) then
                 call la_qgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
                 call la_qlaqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              else
                 call la_qgghd3('N','N',n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,work,- &
                           1,ierr)
                 lwkopt = max(lwkopt,3*n + int(work(1),KIND=ilp))
                 call la_qlaqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alphar,alphai, &
                           beta,vl,ldvl,vr,ldvr,work,-1,0,ierr)
                 lwkopt = max(lwkopt,2*n + int(work(1),KIND=ilp))
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_qlange('M',n,n,a,lda,work)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_qlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_qlange('M',n,n,b,ldb,work)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_qlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           iwrk = iright + n
           call la_qggbal('P',n,a,lda,b,ldb,ilo,ihi,work(ileft),work(iright), &
                     work(iwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = iwrk
           iwrk = itau + irows
           call la_qgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_qormqr('L','T',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_qlaset('FULL',n,n,zero,one,vl,ldvl)
              if (irows > 1) then
                 call la_qlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_qorgqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_qlaset('FULL',n,n,zero,one,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_qgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_qgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_qlaqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alphar,alphai, &
                     beta,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 110
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_qtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 110
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_qggbak('P','L',n,ilo,ihi,work(ileft),work(iright),n,vl, &
                           ldvl,ierr)
                 loop_50: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_50
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vl(jr,jc)) + abs(vl(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_50
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vl(jr,jc) = vl(jr,jc)*temp
                          vl(jr,jc + 1) = vl(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_50
              end if
              if (ilvr) then
                 call la_qggbak('P','R',n,ilo,ihi,work(ileft),work(iright),n,vr, &
                           ldvr,ierr)
                 loop_100: do jc = 1,n
                    if (alphai(jc) < zero) cycle loop_100
                    temp = zero
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)))
                       end do
                    else
                       do jr = 1,n
                          temp = max(temp,abs(vr(jr,jc)) + abs(vr(jr,jc + 1)))
                       end do
                    end if
                    if (temp < smlnum) cycle loop_100
                    temp = one/temp
                    if (alphai(jc) == zero) then
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                       end do
                    else
                       do jr = 1,n
                          vr(jr,jc) = vr(jr,jc)*temp
                          vr(jr,jc + 1) = vr(jr,jc + 1)*temp
                       end do
                    end if
                 end do loop_100
              end if
              ! end of eigenvector calculation
           end if
           ! undo scaling if necessary
           110 continue
           if (ilascl) then
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphar,n,ierr)
              call la_qlascl('G',0,0,anrmto,anrm,n,1,alphai,n,ierr)
           end if
           if (ilbscl) then
              call la_qlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           work(1) = lwkopt
           return
     end subroutine la_qggev3

     !> CGGES: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> CGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_cgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_c) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkmin,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'CGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'CUNMQR',' ',n,1,n,-1))
              if (ilvsl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'CUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -18
           end if
           if (info /= 0) then
              call la_xerbla('CGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_cggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_claset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_cungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_claset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_cgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           call la_chgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: none needed)
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_ctgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_cggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_cggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_clascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_clascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = lwkopt
           return
     end subroutine la_cgges
     !> ZGGES: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> ZGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_zgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_z) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkmin,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'ZGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'ZUNMQR',' ',n,1,n,-1))
              if (ilvsl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'ZUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -18
           end if
           if (info /= 0) then
              call la_xerbla('ZGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_zggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_zlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_zungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_zlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_zgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           call la_zhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: none needed)
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_ztgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_zggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_zggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_zlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_zlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = lwkopt
           return
     end subroutine la_zgges
     !> WGGES: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> WGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_wgges(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_w) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkmin,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'WGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'WUNMQR',' ',n,1,n,-1))
              if (ilvsl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'WUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -18
           end if
           if (info /= 0) then
              call la_xerbla('WGGES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_wggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_wlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_wungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_wlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_wgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           call la_whgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           ! (workspace: none needed)
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_wtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_wggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_wggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_wlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_wlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = lwkopt
           return
     end subroutine la_wgges

     !> CGGESX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the complex Schur form (S,T),
     !> and, optionally, the left and/or right matrices of Schur vectors (VSL
     !> and VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if T is
     !> upper triangular with non-negative diagonal and S is upper
     !> triangular.

     subroutine la_cggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim,alpha, &
      beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,rwork,iwork,liwork,bwork,info)
        use la_constants_sp,only:zero,one,czero,cone

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: rconde(2),rcondv(2),rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_c) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantsb,wantse, &
                     wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,iright,irows, &
                     irwrk,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,smlnum
           ! Local Arrays
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = 2*n
                 maxwrk = n*(1 + la_ilaenv(1,'CGEQRF',' ',n,1,n,0))
                 maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'CUNMQR',' ',n,1,n,-1)))

                 if (ilvsl) then
                    maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'CUNGQR',' ',n,1,n,-1)) &
                              )
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 2
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_cggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_claset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_cungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_claset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_cgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace:    need n)
           iwrk = itau
           call la_chgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 40
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)

              if (ilbscl) call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              ! (complex workspace: if ijob >= 1, need max(1, 2*sdim*(n-sdim))
                                  ! otherwise, need 1 )
              call la_ctgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork,liwork,ierr &
                        )
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -21) then
                  ! not enough complex workspace
                 info = -21
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_cggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_cggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_clascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_clascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           40 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_cggesx
     !> ZGGESX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the complex Schur form (S,T),
     !> and, optionally, the left and/or right matrices of Schur vectors (VSL
     !> and VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if T is
     !> upper triangular with non-negative diagonal and S is upper
     !> triangular.

     subroutine la_zggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim,alpha, &
      beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,rwork,iwork,liwork,bwork,info)
        use la_constants_dp,only:zero,one,czero,cone

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: rconde(2),rcondv(2),rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_z) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantsb,wantse, &
                     wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,iright,irows, &
                     irwrk,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,smlnum
           ! Local Arrays
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = 2*n
                 maxwrk = n*(1 + la_ilaenv(1,'ZGEQRF',' ',n,1,n,0))
                 maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'ZUNMQR',' ',n,1,n,-1)))

                 if (ilvsl) then
                    maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'ZUNGQR',' ',n,1,n,-1)) &
                              )
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 2
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_zggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_zlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_zungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_zlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_zgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace:    need n)
           iwrk = itau
           call la_zhgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 40
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)

              if (ilbscl) call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              ! (complex workspace: if ijob >= 1, need max(1, 2*sdim*(n-sdim))
                                  ! otherwise, need 1 )
              call la_ztgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork,liwork,ierr &
                        )
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -21) then
                  ! not enough complex workspace
                 info = -21
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_zggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_zggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_zlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_zlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           40 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_zggesx
     !> WGGESX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the complex Schur form (S,T),
     !> and, optionally, the left and/or right matrices of Schur vectors (VSL
     !> and VSR).  This gives the generalized Schur factorization
     !> (A,B) = ( (VSL) S (VSR)**H, (VSL) T (VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T; computes
     !> a reciprocal condition number for the average of the selected
     !> eigenvalues (RCONDE); and computes a reciprocal condition number for
     !> the right and left deflating subspaces corresponding to the selected
     !> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
     !> an orthonormal basis for the corresponding left and right eigenspaces
     !> (deflating subspaces).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0 or for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if T is
     !> upper triangular with non-negative diagonal and S is upper
     !> triangular.

     subroutine la_wggesx(jobvsl,jobvsr,sort,selctg,sense,n,a,lda,b,ldb,sdim,alpha, &
      beta,vsl,ldvsl,vsr,ldvsr,rconde,rcondv,work,lwork,rwork,iwork,liwork,bwork,info)
        use la_constants_qp,only:zero,one,czero,cone

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,liwork,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: rconde(2),rcondv(2),rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_w) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantsb,wantse, &
                     wantsn,wantst,wantsv
           integer(ilp) :: i,icols,ierr,ihi,ijob,ijobvl,ijobvr,ileft,ilo,iright,irows, &
                     irwrk,itau,iwrk,liwmin,lwrk,maxwrk,minwrk
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pl,pr,smlnum
           ! Local Arrays
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1 .or. liwork == -1)
           if (wantsn) then
              ijob = 0
           else if (wantse) then
              ijob = 1
           else if (wantsv) then
              ijob = 2
           else if (wantsb) then
              ijob = 4
           end if
           ! test the input arguments
           info = 0
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,n)) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -15
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -17
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv.)
           if (info == 0) then
              if (n > 0) then
                 minwrk = 2*n
                 maxwrk = n*(1 + la_ilaenv(1,'WGEQRF',' ',n,1,n,0))
                 maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'WUNMQR',' ',n,1,n,-1)))

                 if (ilvsl) then
                    maxwrk = max(maxwrk,n*(1 + la_ilaenv(1,'WUNGQR',' ',n,1,n,-1)) &
                              )
                 end if
                 lwrk = maxwrk
                 if (ijob >= 1) lwrk = max(lwrk,n*n/2)
              else
                 minwrk = 1
                 maxwrk = 1
                 lwrk = 1
              end if
              work(1) = lwrk
              if (wantsn .or. n == 0) then
                 liwmin = 1
              else
                 liwmin = n + 2
              end if
              iwork(1) = liwmin
              if (lwork < minwrk .and. .not. lquery) then
                 info = -21
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -24
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGGESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_wggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvsl) then
              call la_wlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_wungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_wlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           call la_wgghrd(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           ! (complex workspace: need n)
           ! (real workspace:    need n)
           iwrk = itau
           call la_whgeqz('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 40
           end if
           ! sort eigenvalues alpha/beta and compute the reciprocal of
           ! condition number(s)
           if (wantst) then
              ! undo scaling on eigenvalues before selctging
              if (ilascl) call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)

              if (ilbscl) call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              ! reorder eigenvalues, transform generalized schur vectors, and
              ! compute reciprocal condition numbers
              ! (complex workspace: if ijob >= 1, need max(1, 2*sdim*(n-sdim))
                                  ! otherwise, need 1 )
              call la_wtgsen(ijob,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pl,pr,dif,work(iwrk),lwork - iwrk + 1,iwork,liwork,ierr &
                        )
              if (ijob >= 1) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (ierr == -21) then
                  ! not enough complex workspace
                 info = -21
              else
                 if (ijob == 1 .or. ijob == 4) then
                    rconde(1) = pl
                    rconde(2) = pr
                 end if
                 if (ijob == 2 .or. ijob == 4) then
                    rcondv(1) = dif(1)
                    rcondv(2) = dif(2)
                 end if
                 if (ierr == 1) info = n + 3
              end if
           end if
           ! apply permutation to vsl and vsr
           ! (workspace: none needed)
           if (ilvsl) call la_wggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_wggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_wlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_wlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           40 continue
           work(1) = maxwrk
           iwork(1) = liwmin
           return
     end subroutine la_wggesx

     !> CGGEV: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_cggev(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkmin,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(sp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
           ! Statement Functions
           real(sp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=sp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'CGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'CUNMQR',' ',n,1,n,0))
              if (ilvl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'CUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -15
           end if
           if (info /= 0) then
              call la_xerbla('CGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('E')*la_slamch('B')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_cggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_claset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_cungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_claset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_cgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_cgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_chgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           ! (real workspace: need 2*n)
           ! (complex workspace: need 2*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_ctgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_cggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_cggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = lwkopt
           return
     end subroutine la_cggev
     !> ZGGEV: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_zggev(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkmin,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(dp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=dp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'ZGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'ZUNMQR',' ',n,1,n,0))
              if (ilvl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'ZUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -15
           end if
           if (info /= 0) then
              call la_xerbla('ZGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('E')*la_dlamch('B')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_zggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_zlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_zungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_zlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_zgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_zgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_zhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           ! (real workspace: need 2*n)
           ! (complex workspace: need 2*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_ztgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_zggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_zggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = lwkopt
           return
     end subroutine la_zggev
     !> WGGEV: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_wggev(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkmin,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(qp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=qp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              lwkmin = max(1,2*n)
              lwkopt = max(1,n + n*la_ilaenv(1,'WGEQRF',' ',n,1,n,0))
              lwkopt = max(lwkopt,n + n*la_ilaenv(1,'WUNMQR',' ',n,1,n,0))
              if (ilvl) then
                 lwkopt = max(lwkopt,n + n*la_ilaenv(1,'WUNGQR',' ',n,1,n,-1))

              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) info = -15
           end if
           if (info /= 0) then
              call la_xerbla('WGGEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('E')*la_qlamch('B')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ! (real workspace: need 6*n)
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_wggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           ! (complex workspace: need n, prefer n*nb)
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           ! (complex workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_wlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_wungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_wlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_wgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_wgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_whgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           ! (real workspace: need 2*n)
           ! (complex workspace: need 2*n)
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_wtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              ! (workspace: none needed)
              if (ilvl) then
                 call la_wggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_wggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = lwkopt
           return
     end subroutine la_wggev

     !> CGGEVX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B) the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> Optionally, it also computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_cggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alpha,beta,vl, &
     ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork,rwork, &
               iwork,bwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(sp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: lscale(*),rconde(*),rcondv(*),rscale(*),rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,wantsb,wantse,wantsn, &
                     wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(sp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
           ! Statement Functions
           real(sp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=sp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (noscl .or. la_lsame(balanc,'S') .or. la_lsame(balanc,'B'))) &
                     then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -13
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -15
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 minwrk = 2*n
                 if (wantse) then
                    minwrk = 4*n
                 else if (wantsv .or. wantsb) then
                    minwrk = 2*n*(n + 1)
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'CGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'CUNMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'CUNGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -25
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (real workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_cggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,rwork,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_clange('1',n,n,a,lda,rwork(1))
           if (ilascl) then
              rwork(1) = abnrm
              call la_slascl('G',0,0,anrmto,anrm,1,1,rwork(1),1,ierr)
              abnrm = rwork(1)
           end if
           bbnrm = la_clange('1',n,n,b,ldb,rwork(1))
           if (ilbscl) then
              rwork(1) = bbnrm
              call la_slascl('G',0,0,bnrmto,bnrm,1,1,rwork(1),1,ierr)
              bbnrm = rwork(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to a
           ! (complex workspace: need n, prefer n*nb)
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_claset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_cungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_claset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_cgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_cgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_chgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 90
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! la_ctgevc: (complex workspace: need 2*n )
                   ! (real workspace:    need 2*n )
           ! la_ctgsna: (complex workspace: need 2*n*n if sense='v' or 'b')
                   ! (integer workspace: need n+2 )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_ctgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work(iwrk),rwork,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 90
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_ctgevc) and estimate condition
                 ! numbers (la_ctgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to
                 ! re-calculate eigenvectors and estimate the condition numbers
                 ! one at a time.
                 do i = 1,n
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    bwork(i) = .true.
                    iwrk = n + 1
                    iwrk1 = iwrk + n
                    if (wantse .or. wantsb) then
                       call la_ctgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,1,m,work(iwrk1),rwork,ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 90
                       end if
                    end if
                    call la_ctgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),1,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                              ierr)
                 end do
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_cggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_50: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vl(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_50
                 temp = one/temp
                 do jr = 1,n
                    vl(jr,jc) = vl(jr,jc)*temp
                 end do
              end do loop_50
           end if
           if (ilvr) then
              call la_cggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_80: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vr(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_80
                 temp = one/temp
                 do jr = 1,n
                    vr(jr,jc) = vr(jr,jc)*temp
                 end do
              end do loop_80
           end if
           ! undo scaling if necessary
           90 continue
           if (ilascl) call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = maxwrk
           return
     end subroutine la_cggevx
     !> ZGGEVX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B) the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> Optionally, it also computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_zggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alpha,beta,vl, &
     ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork,rwork, &
               iwork,bwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(dp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: lscale(*),rconde(*),rcondv(*),rscale(*),rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,wantsb,wantse,wantsn, &
                     wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(dp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=dp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (noscl .or. la_lsame(balanc,'S') .or. la_lsame(balanc,'B'))) &
                     then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -13
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -15
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 minwrk = 2*n
                 if (wantse) then
                    minwrk = 4*n
                 else if (wantsv .or. wantsb) then
                    minwrk = 2*n*(n + 1)
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'ZGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'ZUNMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'ZUNGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -25
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (real workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_zggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,rwork,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_zlange('1',n,n,a,lda,rwork(1))
           if (ilascl) then
              rwork(1) = abnrm
              call la_dlascl('G',0,0,anrmto,anrm,1,1,rwork(1),1,ierr)
              abnrm = rwork(1)
           end if
           bbnrm = la_zlange('1',n,n,b,ldb,rwork(1))
           if (ilbscl) then
              rwork(1) = bbnrm
              call la_dlascl('G',0,0,bnrmto,bnrm,1,1,rwork(1),1,ierr)
              bbnrm = rwork(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to a
           ! (complex workspace: need n, prefer n*nb)
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_zlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_zungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_zlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_zgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_zgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_zhgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 90
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! la_ztgevc: (complex workspace: need 2*n )
                   ! (real workspace:    need 2*n )
           ! la_ztgsna: (complex workspace: need 2*n*n if sense='v' or 'b')
                   ! (integer workspace: need n+2 )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_ztgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work(iwrk),rwork,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 90
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_ztgevc) and estimate condition
                 ! numbers (la_ztgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to
                 ! re-calculate eigenvectors and estimate the condition numbers
                 ! one at a time.
                 do i = 1,n
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    bwork(i) = .true.
                    iwrk = n + 1
                    iwrk1 = iwrk + n
                    if (wantse .or. wantsb) then
                       call la_ztgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,1,m,work(iwrk1),rwork,ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 90
                       end if
                    end if
                    call la_ztgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),1,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                              ierr)
                 end do
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_zggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_50: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vl(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_50
                 temp = one/temp
                 do jr = 1,n
                    vl(jr,jc) = vl(jr,jc)*temp
                 end do
              end do loop_50
           end if
           if (ilvr) then
              call la_zggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_80: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vr(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_80
                 temp = one/temp
                 do jr = 1,n
                    vr(jr,jc) = vr(jr,jc)*temp
                 end do
              end do loop_80
           end if
           ! undo scaling if necessary
           90 continue
           if (ilascl) call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = maxwrk
           return
     end subroutine la_zggevx
     !> WGGEVX: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B) the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> Optionally, it also computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
     !> the eigenvalues (RCONDE), and reciprocal condition numbers for the
     !> right eigenvectors (RCONDV).
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j) .
     !> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
     !> of (A,B) satisfies
     !> u(j)**H * A  = lambda(j) * u(j)**H * B.
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_wggevx(balanc,jobvl,jobvr,sense,n,a,lda,b,ldb,alpha,beta,vl, &
     ldvl,vr,ldvr,ilo,ihi,lscale,rscale,abnrm,bbnrm,rconde,rcondv,work,lwork,rwork, &
               iwork,bwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           real(qp),intent(out) :: abnrm,bbnrm
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: lscale(*),rconde(*),rcondv(*),rscale(*),rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery,noscl,wantsb,wantse,wantsn, &
                     wantsv
           character :: chtemp
           integer(ilp) :: i,icols,ierr,ijobvl,ijobvr,in,irows,itau,iwrk,iwrk1,j,jc, &
                     jr,m,maxwrk,minwrk
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(qp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=qp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           noscl = la_lsame(balanc,'N') .or. la_lsame(balanc,'P')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (.not. (noscl .or. la_lsame(balanc,'S') .or. la_lsame(balanc,'B'))) &
                     then
              info = -1
           else if (ijobvl <= 0) then
              info = -2
           else if (ijobvr <= 0) then
              info = -3
           else if (.not. (wantsn .or. wantse .or. wantsb .or. wantsv)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -13
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -15
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! nb refers to the optimal block size for the immediately
             ! following subroutine, as returned by la_ilaenv. the workspace is
             ! computed assuming ilo = 1 and ihi = n, the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 minwrk = 2*n
                 if (wantse) then
                    minwrk = 4*n
                 else if (wantsv .or. wantsb) then
                    minwrk = 2*n*(n + 1)
                 end if
                 maxwrk = minwrk
                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'WGEQRF',' ',n,1,n,0))

                 maxwrk = max(maxwrk,n + n*la_ilaenv(1,'WUNMQR',' ',n,1,n,0))

                 if (ilvl) then
                    maxwrk = max(maxwrk,n + n*la_ilaenv(1,'WUNGQR',' ',n,1,n,0))

                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -25
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGGEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute and/or balance the matrix pair (a,b)
           ! (real workspace: need 6*n if balanc = 's' or 'b', 1 otherwise)
           call la_wggbal(balanc,n,a,lda,b,ldb,ilo,ihi,lscale,rscale,rwork,ierr)

           ! compute abnrm and bbnrm
           abnrm = la_wlange('1',n,n,a,lda,rwork(1))
           if (ilascl) then
              rwork(1) = abnrm
              call la_qlascl('G',0,0,anrmto,anrm,1,1,rwork(1),1,ierr)
              abnrm = rwork(1)
           end if
           bbnrm = la_wlange('1',n,n,b,ldb,rwork(1))
           if (ilbscl) then
              rwork(1) = bbnrm
              call la_qlascl('G',0,0,bnrmto,bnrm,1,1,rwork(1),1,ierr)
              bbnrm = rwork(1)
           end if
           ! reduce b to triangular form (qr decomposition of b)
           ! (complex workspace: need n, prefer n*nb )
           irows = ihi + 1 - ilo
           if (ilv .or. .not. wantsn) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the unitary transformation to a
           ! (complex workspace: need n, prefer n*nb)
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl and/or vr
           ! (workspace: need n, prefer n*nb)
           if (ilvl) then
              call la_wlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_wungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           if (ilvr) call la_wlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           ! (workspace: none needed)
           if (ilv .or. .not. wantsn) then
              ! eigenvectors requested -- work on whole matrix.
              call la_wgghrd(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        ierr)
           else
              call la_wgghrd('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur forms and schur vectors)
           ! (complex workspace: need n)
           ! (real workspace: need n)
           iwrk = itau
           if (ilv .or. .not. wantsn) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_whgeqz(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 90
           end if
           ! compute eigenvectors and estimate condition numbers if desired
           ! la_wtgevc: (complex workspace: need 2*n )
                   ! (real workspace:    need 2*n )
           ! la_wtgsna: (complex workspace: need 2*n*n if sense='v' or 'b')
                   ! (integer workspace: need n+2 )
           if (ilv .or. .not. wantsn) then
              if (ilv) then
                 if (ilvl) then
                    if (ilvr) then
                       chtemp = 'B'
                    else
                       chtemp = 'L'
                    end if
                 else
                    chtemp = 'R'
                 end if
                 call la_wtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                            in,work(iwrk),rwork,ierr)
                 if (ierr /= 0) then
                    info = n + 2
                    go to 90
                 end if
              end if
              if (.not. wantsn) then
                 ! compute eigenvectors (la_wtgevc) and estimate condition
                 ! numbers (la_wtgsna). note that the definition of the condition
                 ! number is not invariant under transformation (u,v) to
                 ! (q*u, z*v), where (u,v) are eigenvectors of the generalized
                 ! schur form (s,t), q and z are orthogonal matrices. in order
                 ! to avoid using extra 2*n*n workspace, we have to
                 ! re-calculate eigenvectors and estimate the condition numbers
                 ! one at a time.
                 do i = 1,n
                    do j = 1,n
                       bwork(j) = .false.
                    end do
                    bwork(i) = .true.
                    iwrk = n + 1
                    iwrk1 = iwrk + n
                    if (wantse .or. wantsb) then
                       call la_wtgevc('B','S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                                 iwrk),n,1,m,work(iwrk1),rwork,ierr)
                       if (ierr /= 0) then
                          info = n + 2
                          go to 90
                       end if
                    end if
                    call la_wtgsna(sense,'S',bwork,n,a,lda,b,ldb,work(1),n,work( &
                    iwrk),n,rconde(i),rcondv(i),1,m,work(iwrk1),lwork - iwrk1 + 1,iwork, &
                              ierr)
                 end do
              end if
           end if
           ! undo balancing on vl and vr and normalization
           ! (workspace: none needed)
           if (ilvl) then
              call la_wggbak(balanc,'L',n,ilo,ihi,lscale,rscale,n,vl,ldvl,ierr)

              loop_50: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vl(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_50
                 temp = one/temp
                 do jr = 1,n
                    vl(jr,jc) = vl(jr,jc)*temp
                 end do
              end do loop_50
           end if
           if (ilvr) then
              call la_wggbak(balanc,'R',n,ilo,ihi,lscale,rscale,n,vr,ldvr,ierr)

              loop_80: do jc = 1,n
                 temp = zero
                 do jr = 1,n
                    temp = max(temp,abs1(vr(jr,jc)))
                 end do
                 if (temp < smlnum) cycle loop_80
                 temp = one/temp
                 do jr = 1,n
                    vr(jr,jc) = vr(jr,jc)*temp
                 end do
              end do loop_80
           end if
           ! undo scaling if necessary
           90 continue
           if (ilascl) call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = maxwrk
           return
     end subroutine la_wggevx

     !> CGEES: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_cgees(jobvs,sort,select,n,a,lda,sdim,w,vs,ldvs,work,lwork, &
               rwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_c) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,maxwrk, &
                     minwrk
           real(sp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_chseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'CGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_chseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=sp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'CUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_clascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_cgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_cgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_clacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_chseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_clascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (cworkspace: none)
              ! (rworkspace: none)
              call la_ctrsen('N',jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,s,sep,work( &
                        iwrk),lwork - iwrk + 1,icond)
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_cgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_clascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_ccopy(n,a,lda + 1,w,1)
           end if
           work(1) = maxwrk
           return
     end subroutine la_cgees
     !> ZGEES: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_zgees(jobvs,sort,select,n,a,lda,sdim,w,vs,ldvs,work,lwork, &
               rwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_z) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,maxwrk, &
                     minwrk
           real(dp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_zhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'ZGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_zhseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=dp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'ZUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_zlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_zgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_zgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_zlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_zhseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_zlascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (cworkspace: none)
              ! (rworkspace: none)
              call la_ztrsen('N',jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,s,sep,work( &
                        iwrk),lwork - iwrk + 1,icond)
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_zgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_zlascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_zcopy(n,a,lda + 1,w,1)
           end if
           work(1) = maxwrk
           return
     end subroutine la_zgees
     !> WGEES: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left.
     !> The leading columns of Z then form an orthonormal basis for the
     !> invariant subspace corresponding to the selected eigenvalues.
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_wgees(jobvs,sort,select,n,a,lda,sdim,w,vs,ldvs,work,lwork, &
               rwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_w) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantst,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,maxwrk, &
                     minwrk
           real(qp) :: anrm,bignum,cscale,eps,s,sep,smlnum
           ! Local Arrays
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_whseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'WGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_whseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=qp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'WUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGEES ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_wlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_wgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_wgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_wlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_whseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_wlascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues and transform schur vectors
              ! (cworkspace: none)
              ! (rworkspace: none)
              call la_wtrsen('N',jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,s,sep,work( &
                        iwrk),lwork - iwrk + 1,icond)
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_wgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_wlascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_wcopy(n,a,lda + 1,w,1)
           end if
           work(1) = maxwrk
           return
     end subroutine la_wgees

     !> CGEESX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_sp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_cgeesx(jobvs,sort,select,sense,n,a,lda,sdim,w,vs,ldvs,rconde, &
               rcondv,work,lwork,rwork,bwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           real(sp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_c) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantsb,wantse,wantsn,wantst,wantsv,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,lwrk, &
                     maxwrk,minwrk
           real(sp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_chseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_ctrsen later
             ! in the code.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'CGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_chseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=sp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'CUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk, (n*n)/2)
              end if
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_clascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_cgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_cgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_clacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_chseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_clascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (cworkspace: if sense is not 'n', need 2*sdim*(n-sdim)
                           ! otherwise, need none )
              ! (rworkspace: none)
              call la_ctrsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (icond == -14) then
                 ! not enough complex workspace
                 info = -15
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_cgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_clascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_ccopy(n,a,lda + 1,w,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_slascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_cgeesx
     !> ZGEESX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_dp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_zgeesx(jobvs,sort,select,sense,n,a,lda,sdim,w,vs,ldvs,rconde, &
               rcondv,work,lwork,rwork,bwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           real(dp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_z) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantsb,wantse,wantsn,wantst,wantsv,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,lwrk, &
                     maxwrk,minwrk
           real(dp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_zhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_ztrsen later
             ! in the code.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'ZGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_zhseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=dp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'ZUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk, (n*n)/2)
              end if
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_zlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_zgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_zgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_zlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_zhseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_zlascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (cworkspace: if sense is not 'n', need 2*sdim*(n-sdim)
                           ! otherwise, need none )
              ! (rworkspace: none)
              call la_ztrsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (icond == -14) then
                 ! not enough complex workspace
                 info = -15
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_zgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_zlascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_zcopy(n,a,lda + 1,w,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_dlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_zgeesx
     !> WGEESX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
     !> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
     !> Optionally, it also orders the eigenvalues on the diagonal of the
     !> Schur form so that selected eigenvalues are at the top left;
     !> computes a reciprocal condition number for the average of the
     !> selected eigenvalues (RCONDE); and computes a reciprocal condition
     !> number for the right invariant subspace corresponding to the
     !> selected eigenvalues (RCONDV).  The leading columns of Z form an
     !> orthonormal basis for this invariant subspace.
     !> For further explanation of the reciprocal condition numbers RCONDE
     !> and RCONDV, see Section 4.10_qp of the LAPACK Users' Guide (where
     !> these quantities are called s and sep respectively).
     !> A complex matrix is in Schur form if it is upper triangular.

     subroutine la_wgeesx(jobvs,sort,select,sense,n,a,lda,sdim,w,vs,ldvs,rconde, &
               rcondv,work,lwork,rwork,bwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvs,sense,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldvs,lwork,n
           real(qp),intent(out) :: rconde,rcondv
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: vs(ldvs,*),w(*),work(*)
           ! Function Arguments
           procedure(la_select_w) :: select
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantsb,wantse,wantsn,wantst,wantsv,wantvs
           integer(ilp) :: hswork,i,ibal,icond,ierr,ieval,ihi,ilo,itau,iwrk,lwrk, &
                     maxwrk,minwrk
           real(qp) :: anrm,bignum,cscale,eps,smlnum
           ! Local Arrays
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantvs = la_lsame(jobvs,'V')
           wantst = la_lsame(sort,'S')
           wantsn = la_lsame(sense,'N')
           wantse = la_lsame(sense,'E')
           wantsv = la_lsame(sense,'V')
           wantsb = la_lsame(sense,'B')
           lquery = (lwork == -1)
           if ((.not. wantvs) .and. (.not. la_lsame(jobvs,'N'))) then
              info = -1
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -2
           else if (.not. (wantsn .or. wantse .or. wantsv .or. wantsb) .or. (.not. wantst .and. &
                     .not. wantsn)) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
              info = -11
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of real workspace needed at that point in the
             ! code, as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_whseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.
             ! if sense = 'e', 'v' or 'b', then the amount of workspace needed
             ! depends on sdim, which is computed by the routine la_wtrsen later
             ! in the code.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 lwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'WGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 call la_whseqr('S',jobvs,n,1,n,a,lda,w,vs,ldvs,work,-1,ieval)

                 hswork = real(work(1),KIND=qp)
                 if (.not. wantvs) then
                    maxwrk = max(maxwrk,hswork)
                 else
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'WUNGHR',' ',n,1,n,- &
                              1))
                    maxwrk = max(maxwrk,hswork)
                 end if
                 lwrk = maxwrk
                 if (.not. wantsn) lwrk = max(lwrk, (n*n)/2)
              end if
              work(1) = lwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -15
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGEESX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_wlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! permute the matrix to make it more nearly triangular
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_wgebal('P',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = n + itau
           call la_wgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvs) then
              ! copy householder vectors to vs
              call la_wlacpy('L',n,n,a,lda,vs,ldvs)
              ! generate unitary matrix in vs
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vs,ldvs,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
           end if
           sdim = 0
           ! perform qr iteration, accumulating schur vectors in vs if desired
           ! (cworkspace: need 1, prefer hswork (see comments) )
           ! (rworkspace: none)
           iwrk = itau
           call la_whseqr('S',jobvs,n,ilo,ihi,a,lda,w,vs,ldvs,work(iwrk),lwork - &
                     iwrk + 1,ieval)
           if (ieval > 0) info = ieval
           ! sort eigenvalues if desired
           if (wantst .and. info == 0) then
              if (scalea) call la_wlascl('G',0,0,cscale,anrm,n,1,w,n,ierr)
              do i = 1,n
                 bwork(i) = select(w(i))
              end do
              ! reorder eigenvalues, transform schur vectors, and compute
              ! reciprocal condition numbers
              ! (cworkspace: if sense is not 'n', need 2*sdim*(n-sdim)
                           ! otherwise, need none )
              ! (rworkspace: none)
              call la_wtrsen(sense,jobvs,bwork,n,a,lda,vs,ldvs,w,sdim,rconde, &
                        rcondv,work(iwrk),lwork - iwrk + 1,icond)
              if (.not. wantsn) maxwrk = max(maxwrk,2*sdim*(n - sdim))
              if (icond == -14) then
                 ! not enough complex workspace
                 info = -15
              end if
           end if
           if (wantvs) then
              ! undo balancing
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_wgebak('P','R',n,ilo,ihi,rwork(ibal),n,vs,ldvs,ierr)
           end if
           if (scalea) then
              ! undo scaling for the schur form of a
              call la_wlascl('U',0,0,cscale,anrm,n,n,a,lda,ierr)
              call la_wcopy(n,a,lda + 1,w,1)
              if ((wantsv .or. wantsb) .and. info == 0) then
                 dum(1) = rcondv
                 call la_qlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
                 rcondv = dum(1)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_wgeesx

     !> CGEEV: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_cgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork, &
               info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,irwork,itau,iwrk,k,lwork_trevc, &
                     maxwrk,minwrk,nout
           real(sp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(sp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_chseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'CGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 if (wantvl) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'CUNGHR',' ',n,1,n,- &
                              1))
                    call la_ctrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_chseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'CUNGHR',' ',n,1,n,- &
                              1))
                    call la_ctrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_chseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    call la_chseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 end if
                 hswork = int(work(1),KIND=ilp)
                 maxwrk = max(maxwrk,hswork,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_clascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_cgebal('B',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_cgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_clacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_clacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_clacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr('E','N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_chseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need 2*n)
              irwork = ibal + n
              call la_ctrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork(irwork),n,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_cgebak('B','L',n,ilo,ihi,rwork(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_scnrm2(n,vl(1,i),1)
                 call la_csscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vl(k,i),KIND=sp)**2 + aimag(vl(k,i)) &
                              **2
                 end do
                 k = la_isamax(n,rwork(irwork),1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_cscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=sp),zero,KIND=sp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_cgebak('B','R',n,ilo,ihi,rwork(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_scnrm2(n,vr(1,i),1)
                 call la_csscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vr(k,i),KIND=sp)**2 + aimag(vr(k,i)) &
                              **2
                 end do
                 k = la_isamax(n,rwork(irwork),1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_cscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=sp),zero,KIND=sp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_clascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info > 0) then
                 call la_clascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_cgeev
     !> ZGEEV: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_zgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork, &
               info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,irwork,itau,iwrk,k,lwork_trevc, &
                     maxwrk,minwrk,nout
           real(dp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(dp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_zhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'ZGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 if (wantvl) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'ZUNGHR',' ',n,1,n,- &
                              1))
                    call la_ztrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_zhseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'ZUNGHR',' ',n,1,n,- &
                              1))
                    call la_ztrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_zhseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    call la_zhseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 end if
                 hswork = int(work(1),KIND=ilp)
                 maxwrk = max(maxwrk,hswork,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_zlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_zgebal('B',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_zgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_zlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_zlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_zlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr('E','N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_zhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need 2*n)
              irwork = ibal + n
              call la_ztrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork(irwork),n,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_zgebak('B','L',n,ilo,ihi,rwork(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_dznrm2(n,vl(1,i),1)
                 call la_zdscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vl(k,i),KIND=dp)**2 + aimag(vl(k,i)) &
                              **2
                 end do
                 k = la_idamax(n,rwork(irwork),1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_zscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=dp),zero,KIND=dp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_zgebak('B','R',n,ilo,ihi,rwork(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_dznrm2(n,vr(1,i),1)
                 call la_zdscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vr(k,i),KIND=dp)**2 + aimag(vr(k,i)) &
                              **2
                 end do
                 k = la_idamax(n,rwork(irwork),1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_zscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=dp),zero,KIND=dp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_zlascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info > 0) then
                 call la_zlascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_zgeev
     !> WGEEV: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.

     subroutine la_wgeev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork, &
               info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr
           character :: side
           integer(ilp) :: hswork,i,ibal,ierr,ihi,ilo,irwork,itau,iwrk,k,lwork_trevc, &
                     maxwrk,minwrk,nout
           real(qp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(qp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -1
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -10
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_whseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'WGEHRD',' ',n,1,n,0)
                 minwrk = 2*n
                 if (wantvl) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'WUNGHR',' ',n,1,n,- &
                              1))
                    call la_wtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_whseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'WUNGHR',' ',n,1,n,- &
                              1))
                    call la_wtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,n + lwork_trevc)
                    call la_whseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    call la_whseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 end if
                 hswork = int(work(1),KIND=ilp)
                 maxwrk = max(maxwrk,hswork,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGEEV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_wlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix
           ! (cworkspace: none)
           ! (rworkspace: need n)
           ibal = 1
           call la_wgebal('B',n,a,lda,ilo,ihi,rwork(ibal),ierr)
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_wgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_wlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_wlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_wlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr('E','N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_whseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need 2*n)
              irwork = ibal + n
              call la_wtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork(irwork),n,ierr)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_wgebak('B','L',n,ilo,ihi,rwork(ibal),n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_qwnrm2(n,vl(1,i),1)
                 call la_wqscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vl(k,i),KIND=qp)**2 + aimag(vl(k,i)) &
                              **2
                 end do
                 k = la_iqamax(n,rwork(irwork),1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_wscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=qp),zero,KIND=qp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              ! (cworkspace: none)
              ! (rworkspace: need n)
              call la_wgebak('B','R',n,ilo,ihi,rwork(ibal),n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_qwnrm2(n,vr(1,i),1)
                 call la_wqscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(irwork + k - 1) = real(vr(k,i),KIND=qp)**2 + aimag(vr(k,i)) &
                              **2
                 end do
                 k = la_iqamax(n,rwork(irwork),1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(irwork + k - 1))
                 call la_wscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=qp),zero,KIND=qp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_wlascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info > 0) then
                 call la_wlascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_wgeev

     !> CGEEVX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_sp of the LAPACK
     !> Users' Guide.

     subroutine la_cgeevx(balanc,jobvl,jobvr,sense,n,a,lda,w,vl,ldvl,vr,ldvr,ilo, &
               ihi,scale,abnrm,rconde,rcondv,work,lwork,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(sp),intent(out) :: abnrm
           ! Array Arguments
           real(sp),intent(out) :: rconde(*),rcondv(*),rwork(*),scale(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(sp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(sp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') &
                     .or. la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_chseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'CGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_ctrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_chseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    call la_ctrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_chseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    if (wntsnn) then
                       call la_chseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    else
                       call la_chseqr('S','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                 else
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'CUNGHR',' ',n,1,n,- &
                              1))
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,2*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_clange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_clascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_cgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_clange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_slascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_cgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_clacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_clacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_clacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_cunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_chseqr(job,'N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_chseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need n)
              call la_ctrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork,n,ierr)
           end if
           ! compute condition numbers if desired
           ! (cworkspace: need n*n+2*n unless sense = 'e')
           ! (rworkspace: need 2*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_ctrsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,rwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_cgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_scnrm2(n,vl(1,i),1)
                 call la_csscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vl(k,i),KIND=sp)**2 + aimag(vl(k,i))**2
                 end do
                 k = la_isamax(n,rwork,1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(k))
                 call la_cscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=sp),zero,KIND=sp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_cgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_scnrm2(n,vr(1,i),1)
                 call la_csscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vr(k,i),KIND=sp)**2 + aimag(vr(k,i))**2
                 end do
                 k = la_isamax(n,rwork,1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(k))
                 call la_cscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=sp),zero,KIND=sp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_clascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_slascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_clascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_cgeevx
     !> ZGEEVX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_dp of the LAPACK
     !> Users' Guide.

     subroutine la_zgeevx(balanc,jobvl,jobvr,sense,n,a,lda,w,vl,ldvl,vr,ldvr,ilo, &
               ihi,scale,abnrm,rconde,rcondv,work,lwork,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(dp),intent(out) :: abnrm
           ! Array Arguments
           real(dp),intent(out) :: rconde(*),rcondv(*),rwork(*),scale(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(dp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(dp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') &
                     .or. la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_zhseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'ZGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_ztrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_zhseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    call la_ztrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_zhseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    if (wntsnn) then
                       call la_zhseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    else
                       call la_zhseqr('S','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                 else
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'ZUNGHR',' ',n,1,n,- &
                              1))
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,2*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_zlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_zlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_zgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_zlange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_dlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_zgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_zlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_zlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_zlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_zunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_zhseqr(job,'N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_zhseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need n)
              call la_ztrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork,n,ierr)
           end if
           ! compute condition numbers if desired
           ! (cworkspace: need n*n+2*n unless sense = 'e')
           ! (rworkspace: need 2*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_ztrsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,rwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_zgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_dznrm2(n,vl(1,i),1)
                 call la_zdscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vl(k,i),KIND=dp)**2 + aimag(vl(k,i))**2
                 end do
                 k = la_idamax(n,rwork,1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(k))
                 call la_zscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=dp),zero,KIND=dp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_zgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_dznrm2(n,vr(1,i),1)
                 call la_zdscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vr(k,i),KIND=dp)**2 + aimag(vr(k,i))**2
                 end do
                 k = la_idamax(n,rwork,1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(k))
                 call la_zscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=dp),zero,KIND=dp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_zlascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_dlascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_zlascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_zgeevx
     !> WGEEVX: computes for an N-by-N complex nonsymmetric matrix A, the
     !> eigenvalues and, optionally, the left and/or right eigenvectors.
     !> Optionally also, it computes a balancing transformation to improve
     !> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
     !> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
     !> (RCONDE), and reciprocal condition numbers for the right
     !> eigenvectors (RCONDV).
     !> The right eigenvector v(j) of A satisfies
     !> A * v(j) = lambda(j) * v(j)
     !> where lambda(j) is its eigenvalue.
     !> The left eigenvector u(j) of A satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H
     !> where u(j)**H denotes the conjugate transpose of u(j).
     !> The computed eigenvectors are normalized to have Euclidean norm
     !> equal to 1 and largest component real.
     !> Balancing a matrix means permuting the rows and columns to make it
     !> more nearly upper triangular, and applying a diagonal similarity
     !> transformation D * A * D**(-1), where D is a diagonal matrix, to
     !> make its rows and columns closer in norm and the condition numbers
     !> of its eigenvalues and eigenvectors smaller.  The computed
     !> reciprocal condition numbers correspond to the balanced matrix.
     !> Permuting rows and columns will not change the condition numbers
     !> (in exact arithmetic) but diagonal scaling will.  For further
     !> explanation of balancing, see section 4.10.2_qp of the LAPACK
     !> Users' Guide.

     subroutine la_wgeevx(balanc,jobvl,jobvr,sense,n,a,lda,w,vl,ldvl,vr,ldvr,ilo, &
               ihi,scale,abnrm,rconde,rcondv,work,lwork,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: balanc,jobvl,jobvr,sense
           integer(ilp),intent(out) :: ihi,ilo,info
           integer(ilp),intent(in) :: lda,ldvl,ldvr,lwork,n
           real(qp),intent(out) :: abnrm
           ! Array Arguments
           real(qp),intent(out) :: rconde(*),rcondv(*),rwork(*),scale(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: vl(ldvl,*),vr(ldvr,*),w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,scalea,wantvl,wantvr,wntsnb,wntsne,wntsnn,wntsnv
           character :: job,side
           integer(ilp) :: hswork,i,icond,ierr,itau,iwrk,k,lwork_trevc,maxwrk,minwrk, &
                     nout
           real(qp) :: anrm,bignum,cscale,eps,scl,smlnum
           complex(qp) :: tmp
           ! Local Arrays
           logical(lk) :: select(1)
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag,max,sqrt
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           wantvl = la_lsame(jobvl,'V')
           wantvr = la_lsame(jobvr,'V')
           wntsnn = la_lsame(sense,'N')
           wntsne = la_lsame(sense,'E')
           wntsnv = la_lsame(sense,'V')
           wntsnb = la_lsame(sense,'B')
           if (.not. (la_lsame(balanc,'N') .or. la_lsame(balanc,'S') &
                     .or. la_lsame(balanc,'P') .or. la_lsame(balanc,'B'))) then
              info = -1
           else if ((.not. wantvl) .and. (.not. la_lsame(jobvl,'N'))) then
              info = -2
           else if ((.not. wantvr) .and. (.not. la_lsame(jobvr,'N'))) then
              info = -3
           else if (.not. (wntsnn .or. wntsne .or. wntsnb .or. wntsnv) .or. ((wntsne .or. &
                     wntsnb) .and. .not. (wantvl .and. wantvr))) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
              info = -12
           end if
           ! compute workspace
            ! (note: comments in the code beginning "workspace:" describe the
             ! minimal amount of workspace needed at that point in the code,
             ! as well as the preferred amount for good performance.
             ! cworkspace refers to complex workspace, and rworkspace to real
             ! workspace. nb refers to the optimal block size for the
             ! immediately following subroutine, as returned by la_ilaenv.
             ! hswork refers to the workspace preferred by la_whseqr, as
             ! calculated below. hswork is computed assuming ilo=1 and ihi=n,
             ! the worst case.)
           if (info == 0) then
              if (n == 0) then
                 minwrk = 1
                 maxwrk = 1
              else
                 maxwrk = n + n*la_ilaenv(1,'WGEHRD',' ',n,1,n,0)
                 if (wantvl) then
                    call la_wtrevc3('L','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_whseqr('S','V',n,1,n,a,lda,w,vl,ldvl,work,-1,info)

                 else if (wantvr) then
                    call la_wtrevc3('R','B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout, &
                              work,-1,rwork,-1,ierr)
                    lwork_trevc = int(work(1),KIND=ilp)
                    maxwrk = max(maxwrk,lwork_trevc)
                    call la_whseqr('S','V',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                 else
                    if (wntsnn) then
                       call la_whseqr('E','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    else
                       call la_whseqr('S','N',n,1,n,a,lda,w,vr,ldvr,work,-1,info)

                    end if
                 end if
                 hswork = int(work(1),KIND=ilp)
                 if ((.not. wantvl) .and. (.not. wantvr)) then
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                 else
                    minwrk = 2*n
                    if (.not. (wntsnn .or. wntsne)) minwrk = max(minwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,hswork)
                    maxwrk = max(maxwrk,n + (n - 1)*la_ilaenv(1,'WUNGHR',' ',n,1,n,- &
                              1))
                    if (.not. (wntsnn .or. wntsne)) maxwrk = max(maxwrk,n*n + 2*n)
                    maxwrk = max(maxwrk,2*n)
                 end if
                 maxwrk = max(maxwrk,minwrk)
              end if
              work(1) = maxwrk
              if (lwork < minwrk .and. .not. lquery) then
                 info = -20
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGEEVX',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           icond = 0
           anrm = la_wlange('M',n,n,a,lda,dum)
           scalea = .false.
           if (anrm > zero .and. anrm < smlnum) then
              scalea = .true.
              cscale = smlnum
           else if (anrm > bignum) then
              scalea = .true.
              cscale = bignum
           end if
           if (scalea) call la_wlascl('G',0,0,anrm,cscale,n,n,a,lda,ierr)
           ! balance the matrix and compute abnrm
           call la_wgebal(balanc,n,a,lda,ilo,ihi,scale,ierr)
           abnrm = la_wlange('1',n,n,a,lda,dum)
           if (scalea) then
              dum(1) = abnrm
              call la_qlascl('G',0,0,cscale,anrm,1,1,dum,1,ierr)
              abnrm = dum(1)
           end if
           ! reduce to upper hessenberg form
           ! (cworkspace: need 2*n, prefer n+n*nb)
           ! (rworkspace: none)
           itau = 1
           iwrk = itau + n
           call la_wgehrd(n,ilo,ihi,a,lda,work(itau),work(iwrk),lwork - iwrk + 1,ierr &
                     )
           if (wantvl) then
              ! want left eigenvectors
              ! copy householder vectors to vl
              side = 'L'
              call la_wlacpy('L',n,n,a,lda,vl,ldvl)
              ! generate unitary matrix in vl
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vl,ldvl,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vl
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr('S','V',n,ilo,ihi,a,lda,w,vl,ldvl,work(iwrk),lwork - &
                        iwrk + 1,info)
              if (wantvr) then
                 ! want left and right eigenvectors
                 ! copy schur vectors to vr
                 side = 'B'
                 call la_wlacpy('F',n,n,vl,ldvl,vr,ldvr)
              end if
           else if (wantvr) then
              ! want right eigenvectors
              ! copy householder vectors to vr
              side = 'R'
              call la_wlacpy('L',n,n,a,lda,vr,ldvr)
              ! generate unitary matrix in vr
              ! (cworkspace: need 2*n-1, prefer n+(n-1)*nb)
              ! (rworkspace: none)
              call la_wunghr(n,ilo,ihi,vr,ldvr,work(itau),work(iwrk),lwork - iwrk + 1, &
                        ierr)
              ! perform qr iteration, accumulating schur vectors in vr
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr('S','V',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           else
              ! compute eigenvalues only
              ! if condition numbers desired, compute schur form
              if (wntsnn) then
                 job = 'E'
              else
                 job = 'S'
              end if
              ! (cworkspace: need 1, prefer hswork (see comments) )
              ! (rworkspace: none)
              iwrk = itau
              call la_whseqr(job,'N',n,ilo,ihi,a,lda,w,vr,ldvr,work(iwrk),lwork - &
                        iwrk + 1,info)
           end if
           ! if info /= 0 from la_whseqr, then quit
           if (info /= 0) go to 50
           if (wantvl .or. wantvr) then
              ! compute left and/or right eigenvectors
              ! (cworkspace: need 2*n, prefer n + 2*n*nb)
              ! (rworkspace: need n)
              call la_wtrevc3(side,'B',select,n,a,lda,vl,ldvl,vr,ldvr,n,nout,work( &
                         iwrk),lwork - iwrk + 1,rwork,n,ierr)
           end if
           ! compute condition numbers if desired
           ! (cworkspace: need n*n+2*n unless sense = 'e')
           ! (rworkspace: need 2*n unless sense = 'e')
           if (.not. wntsnn) then
              call la_wtrsna(sense,'A',select,n,a,lda,vl,ldvl,vr,ldvr,rconde, &
                        rcondv,n,nout,work(iwrk),n,rwork,icond)
           end if
           if (wantvl) then
              ! undo balancing of left eigenvectors
              call la_wgebak(balanc,'L',n,ilo,ihi,scale,n,vl,ldvl,ierr)
              ! normalize left eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_qwnrm2(n,vl(1,i),1)
                 call la_wqscal(n,scl,vl(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vl(k,i),KIND=qp)**2 + aimag(vl(k,i))**2
                 end do
                 k = la_iqamax(n,rwork,1)
                 tmp = conjg(vl(k,i))/sqrt(rwork(k))
                 call la_wscal(n,tmp,vl(1,i),1)
                 vl(k,i) = cmplx(real(vl(k,i),KIND=qp),zero,KIND=qp)
              end do
           end if
           if (wantvr) then
              ! undo balancing of right eigenvectors
              call la_wgebak(balanc,'R',n,ilo,ihi,scale,n,vr,ldvr,ierr)
              ! normalize right eigenvectors and make largest component real
              do i = 1,n
                 scl = one/la_qwnrm2(n,vr(1,i),1)
                 call la_wqscal(n,scl,vr(1,i),1)
                 do k = 1,n
                    rwork(k) = real(vr(k,i),KIND=qp)**2 + aimag(vr(k,i))**2
                 end do
                 k = la_iqamax(n,rwork,1)
                 tmp = conjg(vr(k,i))/sqrt(rwork(k))
                 call la_wscal(n,tmp,vr(1,i),1)
                 vr(k,i) = cmplx(real(vr(k,i),KIND=qp),zero,KIND=qp)
              end do
           end if
           ! undo scaling if necessary
           50 continue
           if (scalea) then
              call la_wlascl('G',0,0,cscale,anrm,n - info,1,w(info + 1),max(n - info,1) &
                        ,ierr)
              if (info == 0) then
                 if ((wntsnv .or. wntsnb) .and. icond == 0) call la_qlascl('G',0,0,cscale, &
                            anrm,n,1,rcondv,n,ierr)
              else
                 call la_wlascl('G',0,0,cscale,anrm,ilo - 1,1,w,n,ierr)
              end if
           end if
           work(1) = maxwrk
           return
     end subroutine la_wgeevx

     !> CGGES3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> CGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_cgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_c) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(sp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           ! compute workspace
           if (info == 0) then
              call la_cgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,n + int(work(1),KIND=ilp))
              call la_cunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_cungqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              call la_cgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              call la_claqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alpha,beta,vsl, &
                        ldvsl,vsr,ldvsr,work,-1,rwork,0,ierr)
              lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              if (wantst) then
                 call la_ctgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
                           ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)
                 lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('CGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_cggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_claset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_cungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_claset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_cgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_claqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_ctgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_cggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_cggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_clascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_clascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = cmplx(lwkopt,KIND=sp)
           return
     end subroutine la_cgges3
     !> ZGGES3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> ZGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_zgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_z) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(dp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           ! compute workspace
           if (info == 0) then
              call la_zgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,n + int(work(1),KIND=ilp))
              call la_zunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_zungqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              call la_zgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              call la_zlaqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alpha,beta,vsl, &
                        ldvsl,vsr,ldvsr,work,-1,rwork,0,ierr)
              lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              if (wantst) then
                 call la_ztgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
                           ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)
                 lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=dp)
           end if
           if (info /= 0) then
              call la_xerbla('ZGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_zggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_zlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_zungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_zlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_zgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_zlaqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_ztgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_zggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_zggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_zlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_zlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = cmplx(lwkopt,KIND=dp)
           return
     end subroutine la_zgges3
     !> WGGES3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, the generalized complex Schur
     !> form (S, T), and optionally left and/or right Schur vectors (VSL
     !> and VSR). This gives the generalized Schur factorization
     !> (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
     !> where (VSR)**H is the conjugate-transpose of VSR.
     !> Optionally, it also orders the eigenvalues so that a selected cluster
     !> of eigenvalues appears in the leading diagonal blocks of the upper
     !> triangular matrix S and the upper triangular matrix T. The leading
     !> columns of VSL and VSR then form an unitary basis for the
     !> corresponding left and right eigenspaces (deflating subspaces).
     !> (If only the generalized eigenvalues are needed, use the driver
     !> WGGEV instead, which is faster.)
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
     !> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
     !> usually represented as the pair (alpha,beta), as there is a
     !> reasonable interpretation for beta=0, and even for both being zero.
     !> A pair of matrices (S,T) is in generalized complex Schur form if S
     !> and T are upper triangular and, in addition, the diagonal elements
     !> of T are non-negative real numbers.

     subroutine la_wgges3(jobvsl,jobvsr,sort,selctg,n,a,lda,b,ldb,sdim,alpha,beta, &
               vsl,ldvsl,vsr,ldvsr,work,lwork,rwork,bwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvsl,jobvsr,sort
           integer(ilp),intent(out) :: info,sdim
           integer(ilp),intent(in) :: lda,ldb,ldvsl,ldvsr,lwork,n
           ! Array Arguments
           logical(lk),intent(out) :: bwork(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vsl(ldvsl,*),vsr(ldvsr,*),work(*)

           ! Function Arguments
           procedure(la_selctg_w) :: selctg
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: cursl,ilascl,ilbscl,ilvsl,ilvsr,lastsl,lquery,wantst
           integer(ilp) :: i,icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,iright,irows,irwrk, &
                     itau,iwrk,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,pvsl,pvsr,smlnum
           ! Local Arrays
           integer(ilp) :: idum(1)
           real(qp) :: dif(2)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvsl,'N')) then
              ijobvl = 1
              ilvsl = .false.
           else if (la_lsame(jobvsl,'V')) then
              ijobvl = 2
              ilvsl = .true.
           else
              ijobvl = -1
              ilvsl = .false.
           end if
           if (la_lsame(jobvsr,'N')) then
              ijobvr = 1
              ilvsr = .false.
           else if (la_lsame(jobvsr,'V')) then
              ijobvr = 2
              ilvsr = .true.
           else
              ijobvr = -1
              ilvsr = .false.
           end if
           wantst = la_lsame(sort,'S')
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if ((.not. wantst) .and. (.not. la_lsame(sort,'N'))) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldvsl < 1 .or. (ilvsl .and. ldvsl < n)) then
              info = -14
           else if (ldvsr < 1 .or. (ilvsr .and. ldvsr < n)) then
              info = -16
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           ! compute workspace
           if (info == 0) then
              call la_wgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,n + int(work(1),KIND=ilp))
              call la_wunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvsl) then
                 call la_wungqr(n,n,n,vsl,ldvsl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              call la_wgghd3(jobvsl,jobvsr,n,1,n,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                        work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              call la_wlaqz0('S',jobvsl,jobvsr,n,1,n,a,lda,b,ldb,alpha,beta,vsl, &
                        ldvsl,vsr,ldvsr,work,-1,rwork,0,ierr)
              lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              if (wantst) then
                 call la_wtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
                           ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work,-1,idum,1,ierr)
                 lwkopt = max(lwkopt,int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=qp)
           end if
           if (info /= 0) then
              call la_xerbla('WGGES3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) then
              sdim = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrix to make it more nearly triangular
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_wggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           icols = n + 1 - ilo
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vsl
           if (ilvsl) then
              call la_wlaset('FULL',n,n,czero,cone,vsl,ldvsl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vsl(ilo + 1,ilo) &
                           ,ldvsl)
              end if
              call la_wungqr(irows,irows,irows,vsl(ilo,ilo),ldvsl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vsr
           if (ilvsr) call la_wlaset('FULL',n,n,czero,cone,vsr,ldvsr)
           ! reduce to generalized hessenberg form
           call la_wgghd3(jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,vsl,ldvsl,vsr,ldvsr, &
                      work(iwrk),lwork + 1 - iwrk,ierr)
           sdim = 0
           ! perform qz algorithm, computing schur vectors if desired
           iwrk = itau
           call la_wlaqz0('S',jobvsl,jobvsr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vsl, &
                     ldvsl,vsr,ldvsr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 30
           end if
           ! sort eigenvalues alpha/beta if desired
           if (wantst) then
              ! undo scaling on eigenvalues before selecting
              if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,1,alpha,n,ierr)

              if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,1,beta,n,ierr)

              ! select eigenvalues
              do i = 1,n
                 bwork(i) = selctg(alpha(i),beta(i))
              end do
              call la_wtgsen(0,ilvsl,ilvsr,bwork,n,a,lda,b,ldb,alpha,beta,vsl, &
              ldvsl,vsr,ldvsr,sdim,pvsl,pvsr,dif,work(iwrk),lwork - iwrk + 1,idum,1,ierr)

              if (ierr == 1) info = n + 3
           end if
           ! apply back-permutation to vsl and vsr
           if (ilvsl) call la_wggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsl,ldvsl,ierr)
           if (ilvsr) call la_wggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright), &
                     n,vsr,ldvsr,ierr)
           ! undo scaling
           if (ilascl) then
              call la_wlascl('U',0,0,anrmto,anrm,n,n,a,lda,ierr)
              call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           end if
           if (ilbscl) then
              call la_wlascl('U',0,0,bnrmto,bnrm,n,n,b,ldb,ierr)
              call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           end if
           if (wantst) then
              ! check if reordering is correct
              lastsl = .true.
              sdim = 0
              do i = 1,n
                 cursl = selctg(alpha(i),beta(i))
                 if (cursl) sdim = sdim + 1
                 if (cursl .and. .not. lastsl) info = n + 2
                 lastsl = cursl
              end do
           end if
           30 continue
           work(1) = cmplx(lwkopt,KIND=qp)
           return
     end subroutine la_wgges3

     !> CGGEV3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_cggev3(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkopt
           real(sp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(sp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
           ! Statement Functions
           real(sp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=sp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -15
           end if
           ! compute workspace
           if (info == 0) then
              call la_cgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(n,n + int(work(1),KIND=ilp))
              call la_cunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_cungqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              if (ilv) then
                 call la_cgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_claqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              else
                 call la_cgghd3('N','N',n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,work,- &
                           1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_claqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('CGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_slamch('E')*la_slamch('B')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_clange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_clascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_clange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_clascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_cggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_cgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_cunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_claset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_clacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_cungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_claset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_cgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_cgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_claqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_ctgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_cggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_cggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_clascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_clascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = cmplx(lwkopt,KIND=sp)
           return
     end subroutine la_cggev3
     !> ZGGEV3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_zggev3(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkopt
           real(dp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(dp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=dp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -15
           end if
           ! compute workspace
           if (info == 0) then
              call la_zgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,n + int(work(1),KIND=ilp))
              call la_zunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_zungqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              if (ilv) then
                 call la_zgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_zlaqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              else
                 call la_zgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_zlaqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=dp)
           end if
           if (info /= 0) then
              call la_xerbla('ZGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_dlamch('E')*la_dlamch('B')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_zlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_zlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_zlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_zlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_zggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_zgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_zunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_zlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_zlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_zungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_zlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_zgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_zgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_zlaqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_ztgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_zggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_zggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_zlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_zlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = cmplx(lwkopt,KIND=dp)
           return
     end subroutine la_zggev3
     !> WGGEV3: computes for a pair of N-by-N complex nonsymmetric matrices
     !> (A,B), the generalized eigenvalues, and optionally, the left and/or
     !> right generalized eigenvectors.
     !> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
     !> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
     !> singular. It is usually represented as the pair (alpha,beta), as
     !> there is a reasonable interpretation for beta=0, and even for both
     !> being zero.
     !> The right generalized eigenvector v(j) corresponding to the
     !> generalized eigenvalue lambda(j) of (A,B) satisfies
     !> A * v(j) = lambda(j) * B * v(j).
     !> The left generalized eigenvector u(j) corresponding to the
     !> generalized eigenvalues lambda(j) of (A,B) satisfies
     !> u(j)**H * A = lambda(j) * u(j)**H * B
     !> where u(j)**H is the conjugate-transpose of u(j).

     subroutine la_wggev3(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr, &
               work,lwork,rwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobvl,jobvr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldvl,ldvr,lwork,n
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: alpha(*),beta(*),vl(ldvl,*),vr(ldvr,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: ilascl,ilbscl,ilv,ilvl,ilvr,lquery
           character :: chtemp
           integer(ilp) :: icols,ierr,ihi,ijobvl,ijobvr,ileft,ilo,in,iright,irows,irwrk, &
                      itau,iwrk,jc,jr,lwkopt
           real(qp) :: anrm,anrmto,bignum,bnrm,bnrmto,eps,smlnum,temp
           complex(qp) :: x
           ! Local Arrays
           logical(lk) :: ldumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: abs1
           ! Statement Function Definitions
           abs1(x) = abs(real(x,KIND=qp)) + abs(aimag(x))
           ! Executable Statements
           ! decode the input arguments
           if (la_lsame(jobvl,'N')) then
              ijobvl = 1
              ilvl = .false.
           else if (la_lsame(jobvl,'V')) then
              ijobvl = 2
              ilvl = .true.
           else
              ijobvl = -1
              ilvl = .false.
           end if
           if (la_lsame(jobvr,'N')) then
              ijobvr = 1
              ilvr = .false.
           else if (la_lsame(jobvr,'V')) then
              ijobvr = 2
              ilvr = .true.
           else
              ijobvr = -1
              ilvr = .false.
           end if
           ilv = ilvl .or. ilvr
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (ijobvl <= 0) then
              info = -1
           else if (ijobvr <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -15
           end if
           ! compute workspace
           if (info == 0) then
              call la_wgeqrf(n,n,b,ldb,work,work,-1,ierr)
              lwkopt = max(1,n + int(work(1),KIND=ilp))
              call la_wunmqr('L','C',n,n,n,b,ldb,work,a,lda,work,-1,ierr)
              lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              if (ilvl) then
                 call la_wungqr(n,n,n,vl,ldvl,work,work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              if (ilv) then
                 call la_wgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_wlaqz0('S',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              else
                 call la_wgghd3(jobvl,jobvr,n,1,n,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                           work,-1,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
                 call la_wlaqz0('E',jobvl,jobvr,n,1,n,a,lda,b,ldb,alpha,beta,vl, &
                           ldvl,vr,ldvr,work,-1,rwork,0,ierr)
                 lwkopt = max(lwkopt,n + int(work(1),KIND=ilp))
              end if
              work(1) = cmplx(lwkopt,KIND=qp)
           end if
           if (info /= 0) then
              call la_xerbla('WGGEV3 ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! get machine constants
           eps = la_qlamch('E')*la_qlamch('B')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = sqrt(smlnum)/eps
           bignum = one/smlnum
           ! scale a if max element outside range [smlnum,bignum]
           anrm = la_wlange('M',n,n,a,lda,rwork)
           ilascl = .false.
           if (anrm > zero .and. anrm < smlnum) then
              anrmto = smlnum
              ilascl = .true.
           else if (anrm > bignum) then
              anrmto = bignum
              ilascl = .true.
           end if
           if (ilascl) call la_wlascl('G',0,0,anrm,anrmto,n,n,a,lda,ierr)
           ! scale b if max element outside range [smlnum,bignum]
           bnrm = la_wlange('M',n,n,b,ldb,rwork)
           ilbscl = .false.
           if (bnrm > zero .and. bnrm < smlnum) then
              bnrmto = smlnum
              ilbscl = .true.
           else if (bnrm > bignum) then
              bnrmto = bignum
              ilbscl = .true.
           end if
           if (ilbscl) call la_wlascl('G',0,0,bnrm,bnrmto,n,n,b,ldb,ierr)
           ! permute the matrices a, b to isolate eigenvalues if possible
           ileft = 1
           iright = n + 1
           irwrk = iright + n
           call la_wggbal('P',n,a,lda,b,ldb,ilo,ihi,rwork(ileft),rwork(iright), &
                     rwork(irwrk),ierr)
           ! reduce b to triangular form (qr decomposition of b)
           irows = ihi + 1 - ilo
           if (ilv) then
              icols = n + 1 - ilo
           else
              icols = irows
           end if
           itau = 1
           iwrk = itau + irows
           call la_wgeqrf(irows,icols,b(ilo,ilo),ldb,work(itau),work(iwrk),lwork + &
                     1 - iwrk,ierr)
           ! apply the orthogonal transformation to matrix a
           call la_wunmqr('L','C',irows,icols,irows,b(ilo,ilo),ldb,work(itau),a( &
                     ilo,ilo),lda,work(iwrk),lwork + 1 - iwrk,ierr)
           ! initialize vl
           if (ilvl) then
              call la_wlaset('FULL',n,n,czero,cone,vl,ldvl)
              if (irows > 1) then
                 call la_wlacpy('L',irows - 1,irows - 1,b(ilo + 1,ilo),ldb,vl(ilo + 1,ilo), &
                            ldvl)
              end if
              call la_wungqr(irows,irows,irows,vl(ilo,ilo),ldvl,work(itau),work( &
                        iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! initialize vr
           if (ilvr) call la_wlaset('FULL',n,n,czero,cone,vr,ldvr)
           ! reduce to generalized hessenberg form
           if (ilv) then
              ! eigenvectors requested -- work on whole matrix.
              call la_wgghd3(jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,vl,ldvl,vr,ldvr, &
                        work(iwrk),lwork + 1 - iwrk,ierr)
           else
              call la_wgghd3('N','N',irows,1,irows,a(ilo,ilo),lda,b(ilo,ilo), &
                        ldb,vl,ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,ierr)
           end if
           ! perform qz algorithm (compute eigenvalues, and optionally, the
           ! schur form and schur vectors)
           iwrk = itau
           if (ilv) then
              chtemp = 'S'
           else
              chtemp = 'E'
           end if
           call la_wlaqz0(chtemp,jobvl,jobvr,n,ilo,ihi,a,lda,b,ldb,alpha,beta,vl, &
                     ldvl,vr,ldvr,work(iwrk),lwork + 1 - iwrk,rwork(irwrk),0,ierr)
           if (ierr /= 0) then
              if (ierr > 0 .and. ierr <= n) then
                 info = ierr
              else if (ierr > n .and. ierr <= 2*n) then
                 info = ierr - n
              else
                 info = n + 1
              end if
              go to 70
           end if
           ! compute eigenvectors
           if (ilv) then
              if (ilvl) then
                 if (ilvr) then
                    chtemp = 'B'
                 else
                    chtemp = 'L'
                 end if
              else
                 chtemp = 'R'
              end if
              call la_wtgevc(chtemp,'B',ldumma,n,a,lda,b,ldb,vl,ldvl,vr,ldvr,n, &
                        in,work(iwrk),rwork(irwrk),ierr)
              if (ierr /= 0) then
                 info = n + 2
                 go to 70
              end if
              ! undo balancing on vl and vr and normalization
              if (ilvl) then
                 call la_wggbak('P','L',n,ilo,ihi,rwork(ileft),rwork(iright),n,vl, &
                            ldvl,ierr)
                 loop_30: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vl(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_30
                    temp = one/temp
                    do jr = 1,n
                       vl(jr,jc) = vl(jr,jc)*temp
                    end do
                 end do loop_30
              end if
              if (ilvr) then
                 call la_wggbak('P','R',n,ilo,ihi,rwork(ileft),rwork(iright),n,vr, &
                            ldvr,ierr)
                 loop_60: do jc = 1,n
                    temp = zero
                    do jr = 1,n
                       temp = max(temp,abs1(vr(jr,jc)))
                    end do
                    if (temp < smlnum) cycle loop_60
                    temp = one/temp
                    do jr = 1,n
                       vr(jr,jc) = vr(jr,jc)*temp
                    end do
                 end do loop_60
              end if
           end if
           ! undo scaling if necessary
           70 continue
           if (ilascl) call la_wlascl('G',0,0,anrmto,anrm,n,1,alpha,n,ierr)
           if (ilbscl) call la_wlascl('G',0,0,bnrmto,bnrm,n,1,beta,n,ierr)
           work(1) = cmplx(lwkopt,KIND=qp)
           return
     end subroutine la_wggev3

end module la_lapack_eigv_gen
