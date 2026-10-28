!> Nonsymmetric eigenproblem components: Schur factorization, eigenvectors, reordering and condition numbers
module la_lapack_eigv_gen2
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level3_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_eigv_gen3
     use la_lapack_eigv_gen_aux
     use la_lapack_solve_aux
     use la_lapack_solve_tri_comp
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slaein
     public :: la_strevc
     public :: la_strevc3
     public :: la_strsyl
     public :: la_shsein
     public :: la_strsen
     public :: la_strsna
     public :: la_shseqr
     public :: la_dlaein
     public :: la_dtrevc
     public :: la_dtrevc3
     public :: la_dtrsyl
     public :: la_dhsein
     public :: la_dtrsen
     public :: la_dtrsna
     public :: la_dhseqr
#ifdef LA_WITH_XDP
     public :: la_xlaein
     public :: la_xtrevc
     public :: la_xtrevc3
     public :: la_xtrsyl
     public :: la_xhsein
     public :: la_xtrsen
     public :: la_xtrsna
     public :: la_xhseqr
#endif
#ifdef LA_WITH_QP
     public :: la_qlaein
     public :: la_qtrevc
     public :: la_qtrevc3
     public :: la_qtrsyl
     public :: la_qhsein
     public :: la_qtrsen
     public :: la_qtrsna
     public :: la_qhseqr
#endif
     public :: la_ctrevc
     public :: la_ctrevc3
     public :: la_ctrsna
     public :: la_claein
     public :: la_ctrsyl
     public :: la_chsein
     public :: la_ctrsen
     public :: la_chseqr
     public :: la_ztrevc
     public :: la_ztrevc3
     public :: la_ztrsna
     public :: la_zlaein
     public :: la_ztrsyl
     public :: la_zhsein
     public :: la_ztrsen
     public :: la_zhseqr
#ifdef LA_WITH_XDP
     public :: la_ytrevc
     public :: la_ytrevc3
     public :: la_ytrsna
     public :: la_ylaein
     public :: la_ytrsyl
     public :: la_yhsein
     public :: la_ytrsen
     public :: la_yhseqr
#endif
#ifdef LA_WITH_QP
     public :: la_wtrevc
     public :: la_wtrevc3
     public :: la_wtrsna
     public :: la_wlaein
     public :: la_wtrsyl
     public :: la_whsein
     public :: la_wtrsen
     public :: la_whseqr
#endif

     contains

     !> SLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
     !> matrix H.

     pure subroutine la_slaein(rightv,noinit,n,h,ldh,wr,wi,vr,vi,b,ldb,work,eps3, &
               smlnum,bignum,info)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(sp),intent(in) :: bignum,eps3,smlnum,wi,wr
           ! Array Arguments
           real(sp),intent(out) :: b(ldb,*),work(*)
           real(sp),intent(in) :: h(ldh,*)
           real(sp),intent(inout) :: vi(*),vr(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: tenth = 1.0e-1_sp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,i1,i2,i3,ierr,its,j
           real(sp) :: absbii,absbjj,ei,ej,growto,norm,nrmsml,rec,rootn,scale,temp, &
                     vcrit,vmax,vnorm,w,w1,x,xi,xr,y
           ! Intrinsic Functions
           intrinsic :: abs,max,real,sqrt
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=sp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - (wr,wi)*i (except that the subdiagonal elements and
           ! the imaginary parts of the diagonal elements are not stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - wr
           end do
           if (wi == zero) then
              ! real eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                 end do
              else
                 ! scale supplied initial vector.
                 vnorm = la_snrm2(n,vr,1)
                 call la_sscal(n, (eps3*rootn)/max(vnorm,nrmsml),vr,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do i = 1,n - 1
                    ei = h(i + 1,i)
                    if (abs(b(i,i)) < abs(ei)) then
                       ! interchange rows and eliminate.
                       x = b(i,i)/ei
                       b(i,i) = ei
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(i,i) == zero) b(i,i) = eps3
                       x = ei/b(i,i)
                       if (x /= zero) then
                          do j = i + 1,n
                             b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(n,n) == zero) b(n,n) = eps3
                 trans = 'N'
              else
                 ! ul decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do j = n,2,-1
                    ej = h(j,j - 1)
                    if (abs(b(j,j)) < abs(ej)) then
                       ! interchange columns and eliminate.
                       x = b(j,j)/ej
                       b(j,j) = ej
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(j,j) == zero) b(j,j) = eps3
                       x = ej/b(j,j)
                       if (x /= zero) then
                          do i = 1,j - 1
                             b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(1,1) == zero) b(1,1) = eps3
                 trans = 'T'
              end if
              normin = 'N'
              do its = 1,n
                 ! solve u*x = scale*v for a right eigenvector
                   ! or u**t*x = scale*v for a left eigenvector,
                 ! overwriting x on v.
                 call la_slatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,vr,scale,work, &
                            ierr)
                 normin = 'Y'
                 ! test for sufficient growth in the norm of v.
                 vnorm = la_sasum(n,vr,1)
                 if (vnorm >= growto*scale) go to 120
                 ! choose new orthogonal starting vector and try again.
                 temp = eps3/(rootn + one)
                 vr(1) = eps3
                 do i = 2,n
                    vr(i) = temp
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do
              ! failure to find eigenvector in n iterations.
              info = 1
              120 continue
              ! normalize eigenvector.
              i = la_isamax(n,vr,1)
              call la_sscal(n,one/abs(vr(i)),vr,1)
           else
              ! complex eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                    vi(i) = zero
                 end do
              else
                 ! scale supplied initial vector.
                 norm = la_slapy2(la_snrm2(n,vr,1),la_snrm2(n,vi,1))

                 rec = (eps3*rootn)/max(norm,nrmsml)
                 call la_sscal(n,rec,vr,1)
                 call la_sscal(n,rec,vi,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(2,1) = -wi
                 do i = 2,n
                    b(i + 1,1) = zero
                 end do
                 loop_170: do i = 1,n - 1
                    absbii = la_slapy2(b(i,i),b(i + 1,i))
                    ei = h(i + 1,i)
                    if (absbii < abs(ei)) then
                       ! interchange rows and eliminate.
                       xr = b(i,i)/ei
                       xi = b(i + 1,i)/ei
                       b(i,i) = ei
                       b(i + 1,i) = zero
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - xr*temp
                          b(j + 1,i + 1) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(i + 2,i) = -wi
                       b(i + 1,i + 1) = b(i + 1,i + 1) - xi*wi
                       b(i + 2,i + 1) = b(i + 2,i + 1) + xr*wi
                    else
                       ! eliminate without interchanging rows.
                       if (absbii == zero) then
                          b(i,i) = eps3
                          b(i + 1,i) = zero
                          absbii = eps3
                       end if
                       ei = (ei/absbii)/absbii
                       xr = b(i,i)*ei
                       xi = -b(i + 1,i)*ei
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j + 1,i + 1) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(i + 2,i + 1) = b(i + 2,i + 1) - wi
                    end if
                    ! compute 1-norm of offdiagonal elements of i-th row.
                    work(i) = la_sasum(n - i,b(i,i + 1),ldb) + la_sasum(n - i,b(i + 2, &
                              i),1)
                 end do loop_170
                 if (b(n,n) == zero .and. b(n + 1,n) == zero) b(n,n) = eps3
                 work(n) = zero
                 i1 = n
                 i2 = 1
                 i3 = -1
              else
                 ! ul decomposition with partial pivoting of conjg(b),
                 ! replacing zero pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(n + 1,n) = wi
                 do j = 1,n - 1
                    b(n + 1,j) = zero
                 end do
                 loop_210: do j = n,2,-1
                    ej = h(j,j - 1)
                    absbjj = la_slapy2(b(j,j),b(j + 1,j))
                    if (absbjj < abs(ej)) then
                       ! interchange columns and eliminate
                       xr = b(j,j)/ej
                       xi = b(j + 1,j)/ej
                       b(j,j) = ej
                       b(j + 1,j) = zero
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - xr*temp
                          b(j,i) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(j + 1,j - 1) = wi
                       b(j - 1,j - 1) = b(j - 1,j - 1) + xi*wi
                       b(j,j - 1) = b(j,j - 1) - xr*wi
                    else
                       ! eliminate without interchange.
                       if (absbjj == zero) then
                          b(j,j) = eps3
                          b(j + 1,j) = zero
                          absbjj = eps3
                       end if
                       ej = (ej/absbjj)/absbjj
                       xr = b(j,j)*ej
                       xi = -b(j + 1,j)*ej
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j,i) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(j,j - 1) = b(j,j - 1) + wi
                    end if
                    ! compute 1-norm of offdiagonal elements of j-th column.
                    work(j) = la_sasum(j - 1,b(1,j),1) + la_sasum(j - 1,b(j + 1,1), &
                               ldb)
                 end do loop_210
                 if (b(1,1) == zero .and. b(2,1) == zero) b(1,1) = eps3
                 work(1) = zero
                 i1 = 1
                 i2 = n
                 i3 = 1
              end if
              loop_270: do its = 1,n
                 scale = one
                 vmax = one
                 vcrit = bignum
                 ! solve u*(xr,xi) = scale*(vr,vi) for a right eigenvector,
                   ! or u**t*(xr,xi) = scale*(vr,vi) for a left eigenvector,
                 ! overwriting (xr,xi) on (vr,vi).
                 loop_250: do i = i1,i2,i3
                    if (work(i) > vcrit) then
                       rec = one/vmax
                       call la_sscal(n,rec,vr,1)
                       call la_sscal(n,rec,vi,1)
                       scale = scale*rec
                       vmax = one
                       vcrit = bignum
                    end if
                    xr = vr(i)
                    xi = vi(i)
                    if (rightv) then
                       do j = i + 1,n
                          xr = xr - b(i,j)*vr(j) + b(j + 1,i)*vi(j)
                          xi = xi - b(i,j)*vi(j) - b(j + 1,i)*vr(j)
                       end do
                    else
                       do j = 1,i - 1
                          xr = xr - b(j,i)*vr(j) + b(i + 1,j)*vi(j)
                          xi = xi - b(j,i)*vi(j) - b(i + 1,j)*vr(j)
                       end do
                    end if
                    w = abs(b(i,i)) + abs(b(i + 1,i))
                    if (w > smlnum) then
                       if (w < one) then
                          w1 = abs(xr) + abs(xi)
                          if (w1 > w*bignum) then
                             rec = one/w1
                             call la_sscal(n,rec,vr,1)
                             call la_sscal(n,rec,vi,1)
                             xr = vr(i)
                             xi = vi(i)
                             scale = scale*rec
                             vmax = vmax*rec
                          end if
                       end if
                       ! divide by diagonal element of b.
                       call la_sladiv(xr,xi,b(i,i),b(i + 1,i),vr(i),vi(i))

                       vmax = max(abs(vr(i)) + abs(vi(i)),vmax)
                       vcrit = bignum/vmax
                    else
                       do j = 1,n
                          vr(j) = zero
                          vi(j) = zero
                       end do
                       vr(i) = one
                       vi(i) = one
                       scale = zero
                       vmax = one
                       vcrit = bignum
                    end if
                 end do loop_250
                 ! test for sufficient growth in the norm of (vr,vi).
                 vnorm = la_sasum(n,vr,1) + la_sasum(n,vi,1)
                 if (vnorm >= growto*scale) go to 280
                 ! choose a new orthogonal starting vector and try again.
                 y = eps3/(rootn + one)
                 vr(1) = eps3
                 vi(1) = zero
                 do i = 2,n
                    vr(i) = y
                    vi(i) = zero
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do loop_270
              ! failure to find eigenvector in n iterations
              info = 1
              280 continue
              ! normalize eigenvector.
              vnorm = zero
              do i = 1,n
                 vnorm = max(vnorm,abs(vr(i)) + abs(vi(i)))
              end do
              call la_sscal(n,one/vnorm,vr,1)
              call la_sscal(n,one/vnorm,vi,1)
           end if
           return
     end subroutine la_slaein
     !> DLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
     !> matrix H.

     pure subroutine la_dlaein(rightv,noinit,n,h,ldh,wr,wi,vr,vi,b,ldb,work,eps3, &
               smlnum,bignum,info)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(dp),intent(in) :: bignum,eps3,smlnum,wi,wr
           ! Array Arguments
           real(dp),intent(out) :: b(ldb,*),work(*)
           real(dp),intent(in) :: h(ldh,*)
           real(dp),intent(inout) :: vi(*),vr(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: tenth = 1.0e-1_dp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,i1,i2,i3,ierr,its,j
           real(dp) :: absbii,absbjj,ei,ej,growto,norm,nrmsml,rec,rootn,scale,temp, &
                     vcrit,vmax,vnorm,w,w1,x,xi,xr,y
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=dp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - (wr,wi)*i (except that the subdiagonal elements and
           ! the imaginary parts of the diagonal elements are not stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - wr
           end do
           if (wi == zero) then
              ! real eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                 end do
              else
                 ! scale supplied initial vector.
                 vnorm = la_dnrm2(n,vr,1)
                 call la_dscal(n, (eps3*rootn)/max(vnorm,nrmsml),vr,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do i = 1,n - 1
                    ei = h(i + 1,i)
                    if (abs(b(i,i)) < abs(ei)) then
                       ! interchange rows and eliminate.
                       x = b(i,i)/ei
                       b(i,i) = ei
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(i,i) == zero) b(i,i) = eps3
                       x = ei/b(i,i)
                       if (x /= zero) then
                          do j = i + 1,n
                             b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(n,n) == zero) b(n,n) = eps3
                 trans = 'N'
              else
                 ! ul decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do j = n,2,-1
                    ej = h(j,j - 1)
                    if (abs(b(j,j)) < abs(ej)) then
                       ! interchange columns and eliminate.
                       x = b(j,j)/ej
                       b(j,j) = ej
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(j,j) == zero) b(j,j) = eps3
                       x = ej/b(j,j)
                       if (x /= zero) then
                          do i = 1,j - 1
                             b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(1,1) == zero) b(1,1) = eps3
                 trans = 'T'
              end if
              normin = 'N'
              do its = 1,n
                 ! solve u*x = scale*v for a right eigenvector
                   ! or u**t*x = scale*v for a left eigenvector,
                 ! overwriting x on v.
                 call la_dlatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,vr,scale,work, &
                            ierr)
                 normin = 'Y'
                 ! test for sufficient growth in the norm of v.
                 vnorm = la_dasum(n,vr,1)
                 if (vnorm >= growto*scale) go to 120
                 ! choose new orthogonal starting vector and try again.
                 temp = eps3/(rootn + one)
                 vr(1) = eps3
                 do i = 2,n
                    vr(i) = temp
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do
              ! failure to find eigenvector in n iterations.
              info = 1
              120 continue
              ! normalize eigenvector.
              i = la_idamax(n,vr,1)
              call la_dscal(n,one/abs(vr(i)),vr,1)
           else
              ! complex eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                    vi(i) = zero
                 end do
              else
                 ! scale supplied initial vector.
                 norm = la_dlapy2(la_dnrm2(n,vr,1),la_dnrm2(n,vi,1))

                 rec = (eps3*rootn)/max(norm,nrmsml)
                 call la_dscal(n,rec,vr,1)
                 call la_dscal(n,rec,vi,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(2,1) = -wi
                 do i = 2,n
                    b(i + 1,1) = zero
                 end do
                 loop_170: do i = 1,n - 1
                    absbii = la_dlapy2(b(i,i),b(i + 1,i))
                    ei = h(i + 1,i)
                    if (absbii < abs(ei)) then
                       ! interchange rows and eliminate.
                       xr = b(i,i)/ei
                       xi = b(i + 1,i)/ei
                       b(i,i) = ei
                       b(i + 1,i) = zero
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - xr*temp
                          b(j + 1,i + 1) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(i + 2,i) = -wi
                       b(i + 1,i + 1) = b(i + 1,i + 1) - xi*wi
                       b(i + 2,i + 1) = b(i + 2,i + 1) + xr*wi
                    else
                       ! eliminate without interchanging rows.
                       if (absbii == zero) then
                          b(i,i) = eps3
                          b(i + 1,i) = zero
                          absbii = eps3
                       end if
                       ei = (ei/absbii)/absbii
                       xr = b(i,i)*ei
                       xi = -b(i + 1,i)*ei
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j + 1,i + 1) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(i + 2,i + 1) = b(i + 2,i + 1) - wi
                    end if
                    ! compute 1-norm of offdiagonal elements of i-th row.
                    work(i) = la_dasum(n - i,b(i,i + 1),ldb) + la_dasum(n - i,b(i + 2, &
                              i),1)
                 end do loop_170
                 if (b(n,n) == zero .and. b(n + 1,n) == zero) b(n,n) = eps3
                 work(n) = zero
                 i1 = n
                 i2 = 1
                 i3 = -1
              else
                 ! ul decomposition with partial pivoting of conjg(b),
                 ! replacing zero pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(n + 1,n) = wi
                 do j = 1,n - 1
                    b(n + 1,j) = zero
                 end do
                 loop_210: do j = n,2,-1
                    ej = h(j,j - 1)
                    absbjj = la_dlapy2(b(j,j),b(j + 1,j))
                    if (absbjj < abs(ej)) then
                       ! interchange columns and eliminate
                       xr = b(j,j)/ej
                       xi = b(j + 1,j)/ej
                       b(j,j) = ej
                       b(j + 1,j) = zero
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - xr*temp
                          b(j,i) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(j + 1,j - 1) = wi
                       b(j - 1,j - 1) = b(j - 1,j - 1) + xi*wi
                       b(j,j - 1) = b(j,j - 1) - xr*wi
                    else
                       ! eliminate without interchange.
                       if (absbjj == zero) then
                          b(j,j) = eps3
                          b(j + 1,j) = zero
                          absbjj = eps3
                       end if
                       ej = (ej/absbjj)/absbjj
                       xr = b(j,j)*ej
                       xi = -b(j + 1,j)*ej
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j,i) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(j,j - 1) = b(j,j - 1) + wi
                    end if
                    ! compute 1-norm of offdiagonal elements of j-th column.
                    work(j) = la_dasum(j - 1,b(1,j),1) + la_dasum(j - 1,b(j + 1,1), &
                               ldb)
                 end do loop_210
                 if (b(1,1) == zero .and. b(2,1) == zero) b(1,1) = eps3
                 work(1) = zero
                 i1 = 1
                 i2 = n
                 i3 = 1
              end if
              loop_270: do its = 1,n
                 scale = one
                 vmax = one
                 vcrit = bignum
                 ! solve u*(xr,xi) = scale*(vr,vi) for a right eigenvector,
                   ! or u**t*(xr,xi) = scale*(vr,vi) for a left eigenvector,
                 ! overwriting (xr,xi) on (vr,vi).
                 loop_250: do i = i1,i2,i3
                    if (work(i) > vcrit) then
                       rec = one/vmax
                       call la_dscal(n,rec,vr,1)
                       call la_dscal(n,rec,vi,1)
                       scale = scale*rec
                       vmax = one
                       vcrit = bignum
                    end if
                    xr = vr(i)
                    xi = vi(i)
                    if (rightv) then
                       do j = i + 1,n
                          xr = xr - b(i,j)*vr(j) + b(j + 1,i)*vi(j)
                          xi = xi - b(i,j)*vi(j) - b(j + 1,i)*vr(j)
                       end do
                    else
                       do j = 1,i - 1
                          xr = xr - b(j,i)*vr(j) + b(i + 1,j)*vi(j)
                          xi = xi - b(j,i)*vi(j) - b(i + 1,j)*vr(j)
                       end do
                    end if
                    w = abs(b(i,i)) + abs(b(i + 1,i))
                    if (w > smlnum) then
                       if (w < one) then
                          w1 = abs(xr) + abs(xi)
                          if (w1 > w*bignum) then
                             rec = one/w1
                             call la_dscal(n,rec,vr,1)
                             call la_dscal(n,rec,vi,1)
                             xr = vr(i)
                             xi = vi(i)
                             scale = scale*rec
                             vmax = vmax*rec
                          end if
                       end if
                       ! divide by diagonal element of b.
                       call la_dladiv(xr,xi,b(i,i),b(i + 1,i),vr(i),vi(i))

                       vmax = max(abs(vr(i)) + abs(vi(i)),vmax)
                       vcrit = bignum/vmax
                    else
                       do j = 1,n
                          vr(j) = zero
                          vi(j) = zero
                       end do
                       vr(i) = one
                       vi(i) = one
                       scale = zero
                       vmax = one
                       vcrit = bignum
                    end if
                 end do loop_250
                 ! test for sufficient growth in the norm of (vr,vi).
                 vnorm = la_dasum(n,vr,1) + la_dasum(n,vi,1)
                 if (vnorm >= growto*scale) go to 280
                 ! choose a new orthogonal starting vector and try again.
                 y = eps3/(rootn + one)
                 vr(1) = eps3
                 vi(1) = zero
                 do i = 2,n
                    vr(i) = y
                    vi(i) = zero
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do loop_270
              ! failure to find eigenvector in n iterations
              info = 1
              280 continue
              ! normalize eigenvector.
              vnorm = zero
              do i = 1,n
                 vnorm = max(vnorm,abs(vr(i)) + abs(vi(i)))
              end do
              call la_dscal(n,one/vnorm,vr,1)
              call la_dscal(n,one/vnorm,vi,1)
           end if
           return
     end subroutine la_dlaein
#ifdef LA_WITH_XDP
     !> XLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
     !> matrix H.

     pure subroutine la_xlaein(rightv,noinit,n,h,ldh,wr,wi,vr,vi,b,ldb,work,eps3, &
               smlnum,bignum,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(xdp),intent(in) :: bignum,eps3,smlnum,wi,wr
           ! Array Arguments
           real(xdp),intent(out) :: b(ldb,*),work(*)
           real(xdp),intent(in) :: h(ldh,*)
           real(xdp),intent(inout) :: vi(*),vr(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: tenth = 1.0e-1_xdp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,i1,i2,i3,ierr,its,j
           real(xdp) :: absbii,absbjj,ei,ej,growto,norm,nrmsml,rec,rootn,scale,temp, &
                     vcrit,vmax,vnorm,w,w1,x,xi,xr,y
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=xdp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - (wr,wi)*i (except that the subdiagonal elements and
           ! the imaginary parts of the diagonal elements are not stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - wr
           end do
           if (wi == zero) then
              ! real eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                 end do
              else
                 ! scale supplied initial vector.
                 vnorm = la_xnrm2(n,vr,1)
                 call la_xscal(n, (eps3*rootn)/max(vnorm,nrmsml),vr,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do i = 1,n - 1
                    ei = h(i + 1,i)
                    if (abs(b(i,i)) < abs(ei)) then
                       ! interchange rows and eliminate.
                       x = b(i,i)/ei
                       b(i,i) = ei
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(i,i) == zero) b(i,i) = eps3
                       x = ei/b(i,i)
                       if (x /= zero) then
                          do j = i + 1,n
                             b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(n,n) == zero) b(n,n) = eps3
                 trans = 'N'
              else
                 ! ul decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do j = n,2,-1
                    ej = h(j,j - 1)
                    if (abs(b(j,j)) < abs(ej)) then
                       ! interchange columns and eliminate.
                       x = b(j,j)/ej
                       b(j,j) = ej
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(j,j) == zero) b(j,j) = eps3
                       x = ej/b(j,j)
                       if (x /= zero) then
                          do i = 1,j - 1
                             b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(1,1) == zero) b(1,1) = eps3
                 trans = 'T'
              end if
              normin = 'N'
              do its = 1,n
                 ! solve u*x = scale*v for a right eigenvector
                   ! or u**t*x = scale*v for a left eigenvector,
                 ! overwriting x on v.
                 call la_xlatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,vr,scale,work, &
                            ierr)
                 normin = 'Y'
                 ! test for sufficient growth in the norm of v.
                 vnorm = la_xasum(n,vr,1)
                 if (vnorm >= growto*scale) go to 120
                 ! choose new orthogonal starting vector and try again.
                 temp = eps3/(rootn + one)
                 vr(1) = eps3
                 do i = 2,n
                    vr(i) = temp
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do
              ! failure to find eigenvector in n iterations.
              info = 1
              120 continue
              ! normalize eigenvector.
              i = la_ixamax(n,vr,1)
              call la_xscal(n,one/abs(vr(i)),vr,1)
           else
              ! complex eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                    vi(i) = zero
                 end do
              else
                 ! scale supplied initial vector.
                 norm = la_xlapy2(la_xnrm2(n,vr,1),la_xnrm2(n,vi,1))

                 rec = (eps3*rootn)/max(norm,nrmsml)
                 call la_xscal(n,rec,vr,1)
                 call la_xscal(n,rec,vi,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(2,1) = -wi
                 do i = 2,n
                    b(i + 1,1) = zero
                 end do
                 loop_170: do i = 1,n - 1
                    absbii = la_xlapy2(b(i,i),b(i + 1,i))
                    ei = h(i + 1,i)
                    if (absbii < abs(ei)) then
                       ! interchange rows and eliminate.
                       xr = b(i,i)/ei
                       xi = b(i + 1,i)/ei
                       b(i,i) = ei
                       b(i + 1,i) = zero
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - xr*temp
                          b(j + 1,i + 1) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(i + 2,i) = -wi
                       b(i + 1,i + 1) = b(i + 1,i + 1) - xi*wi
                       b(i + 2,i + 1) = b(i + 2,i + 1) + xr*wi
                    else
                       ! eliminate without interchanging rows.
                       if (absbii == zero) then
                          b(i,i) = eps3
                          b(i + 1,i) = zero
                          absbii = eps3
                       end if
                       ei = (ei/absbii)/absbii
                       xr = b(i,i)*ei
                       xi = -b(i + 1,i)*ei
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j + 1,i + 1) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(i + 2,i + 1) = b(i + 2,i + 1) - wi
                    end if
                    ! compute 1-norm of offdiagonal elements of i-th row.
                    work(i) = la_xasum(n - i,b(i,i + 1),ldb) + la_xasum(n - i,b(i + 2, &
                              i),1)
                 end do loop_170
                 if (b(n,n) == zero .and. b(n + 1,n) == zero) b(n,n) = eps3
                 work(n) = zero
                 i1 = n
                 i2 = 1
                 i3 = -1
              else
                 ! ul decomposition with partial pivoting of conjg(b),
                 ! replacing zero pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(n + 1,n) = wi
                 do j = 1,n - 1
                    b(n + 1,j) = zero
                 end do
                 loop_210: do j = n,2,-1
                    ej = h(j,j - 1)
                    absbjj = la_xlapy2(b(j,j),b(j + 1,j))
                    if (absbjj < abs(ej)) then
                       ! interchange columns and eliminate
                       xr = b(j,j)/ej
                       xi = b(j + 1,j)/ej
                       b(j,j) = ej
                       b(j + 1,j) = zero
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - xr*temp
                          b(j,i) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(j + 1,j - 1) = wi
                       b(j - 1,j - 1) = b(j - 1,j - 1) + xi*wi
                       b(j,j - 1) = b(j,j - 1) - xr*wi
                    else
                       ! eliminate without interchange.
                       if (absbjj == zero) then
                          b(j,j) = eps3
                          b(j + 1,j) = zero
                          absbjj = eps3
                       end if
                       ej = (ej/absbjj)/absbjj
                       xr = b(j,j)*ej
                       xi = -b(j + 1,j)*ej
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j,i) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(j,j - 1) = b(j,j - 1) + wi
                    end if
                    ! compute 1-norm of offdiagonal elements of j-th column.
                    work(j) = la_xasum(j - 1,b(1,j),1) + la_xasum(j - 1,b(j + 1,1), &
                               ldb)
                 end do loop_210
                 if (b(1,1) == zero .and. b(2,1) == zero) b(1,1) = eps3
                 work(1) = zero
                 i1 = 1
                 i2 = n
                 i3 = 1
              end if
              loop_270: do its = 1,n
                 scale = one
                 vmax = one
                 vcrit = bignum
                 ! solve u*(xr,xi) = scale*(vr,vi) for a right eigenvector,
                   ! or u**t*(xr,xi) = scale*(vr,vi) for a left eigenvector,
                 ! overwriting (xr,xi) on (vr,vi).
                 loop_250: do i = i1,i2,i3
                    if (work(i) > vcrit) then
                       rec = one/vmax
                       call la_xscal(n,rec,vr,1)
                       call la_xscal(n,rec,vi,1)
                       scale = scale*rec
                       vmax = one
                       vcrit = bignum
                    end if
                    xr = vr(i)
                    xi = vi(i)
                    if (rightv) then
                       do j = i + 1,n
                          xr = xr - b(i,j)*vr(j) + b(j + 1,i)*vi(j)
                          xi = xi - b(i,j)*vi(j) - b(j + 1,i)*vr(j)
                       end do
                    else
                       do j = 1,i - 1
                          xr = xr - b(j,i)*vr(j) + b(i + 1,j)*vi(j)
                          xi = xi - b(j,i)*vi(j) - b(i + 1,j)*vr(j)
                       end do
                    end if
                    w = abs(b(i,i)) + abs(b(i + 1,i))
                    if (w > smlnum) then
                       if (w < one) then
                          w1 = abs(xr) + abs(xi)
                          if (w1 > w*bignum) then
                             rec = one/w1
                             call la_xscal(n,rec,vr,1)
                             call la_xscal(n,rec,vi,1)
                             xr = vr(i)
                             xi = vi(i)
                             scale = scale*rec
                             vmax = vmax*rec
                          end if
                       end if
                       ! divide by diagonal element of b.
                       call la_xladiv(xr,xi,b(i,i),b(i + 1,i),vr(i),vi(i))

                       vmax = max(abs(vr(i)) + abs(vi(i)),vmax)
                       vcrit = bignum/vmax
                    else
                       do j = 1,n
                          vr(j) = zero
                          vi(j) = zero
                       end do
                       vr(i) = one
                       vi(i) = one
                       scale = zero
                       vmax = one
                       vcrit = bignum
                    end if
                 end do loop_250
                 ! test for sufficient growth in the norm of (vr,vi).
                 vnorm = la_xasum(n,vr,1) + la_xasum(n,vi,1)
                 if (vnorm >= growto*scale) go to 280
                 ! choose a new orthogonal starting vector and try again.
                 y = eps3/(rootn + one)
                 vr(1) = eps3
                 vi(1) = zero
                 do i = 2,n
                    vr(i) = y
                    vi(i) = zero
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do loop_270
              ! failure to find eigenvector in n iterations
              info = 1
              280 continue
              ! normalize eigenvector.
              vnorm = zero
              do i = 1,n
                 vnorm = max(vnorm,abs(vr(i)) + abs(vi(i)))
              end do
              call la_xscal(n,one/vnorm,vr,1)
              call la_xscal(n,one/vnorm,vi,1)
           end if
           return
     end subroutine la_xlaein
#endif
#ifdef LA_WITH_QP
     !> QLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
     !> matrix H.

     pure subroutine la_qlaein(rightv,noinit,n,h,ldh,wr,wi,vr,vi,b,ldb,work,eps3, &
               smlnum,bignum,info)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(qp),intent(in) :: bignum,eps3,smlnum,wi,wr
           ! Array Arguments
           real(qp),intent(out) :: b(ldb,*),work(*)
           real(qp),intent(in) :: h(ldh,*)
           real(qp),intent(inout) :: vi(*),vr(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: tenth = 1.0e-1_qp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,i1,i2,i3,ierr,its,j
           real(qp) :: absbii,absbjj,ei,ej,growto,norm,nrmsml,rec,rootn,scale,temp, &
                     vcrit,vmax,vnorm,w,w1,x,xi,xr,y
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=qp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - (wr,wi)*i (except that the subdiagonal elements and
           ! the imaginary parts of the diagonal elements are not stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - wr
           end do
           if (wi == zero) then
              ! real eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                 end do
              else
                 ! scale supplied initial vector.
                 vnorm = la_qnrm2(n,vr,1)
                 call la_qscal(n, (eps3*rootn)/max(vnorm,nrmsml),vr,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do i = 1,n - 1
                    ei = h(i + 1,i)
                    if (abs(b(i,i)) < abs(ei)) then
                       ! interchange rows and eliminate.
                       x = b(i,i)/ei
                       b(i,i) = ei
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(i,i) == zero) b(i,i) = eps3
                       x = ei/b(i,i)
                       if (x /= zero) then
                          do j = i + 1,n
                             b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(n,n) == zero) b(n,n) = eps3
                 trans = 'N'
              else
                 ! ul decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 do j = n,2,-1
                    ej = h(j,j - 1)
                    if (abs(b(j,j)) < abs(ej)) then
                       ! interchange columns and eliminate.
                       x = b(j,j)/ej
                       b(j,j) = ej
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - x*temp
                          b(i,j) = temp
                       end do
                    else
                       ! eliminate without interchange.
                       if (b(j,j) == zero) b(j,j) = eps3
                       x = ej/b(j,j)
                       if (x /= zero) then
                          do i = 1,j - 1
                             b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                          end do
                       end if
                    end if
                 end do
                 if (b(1,1) == zero) b(1,1) = eps3
                 trans = 'T'
              end if
              normin = 'N'
              do its = 1,n
                 ! solve u*x = scale*v for a right eigenvector
                   ! or u**t*x = scale*v for a left eigenvector,
                 ! overwriting x on v.
                 call la_qlatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,vr,scale,work, &
                            ierr)
                 normin = 'Y'
                 ! test for sufficient growth in the norm of v.
                 vnorm = la_qasum(n,vr,1)
                 if (vnorm >= growto*scale) go to 120
                 ! choose new orthogonal starting vector and try again.
                 temp = eps3/(rootn + one)
                 vr(1) = eps3
                 do i = 2,n
                    vr(i) = temp
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do
              ! failure to find eigenvector in n iterations.
              info = 1
              120 continue
              ! normalize eigenvector.
              i = la_iqamax(n,vr,1)
              call la_qscal(n,one/abs(vr(i)),vr,1)
           else
              ! complex eigenvalue.
              if (noinit) then
                 ! set initial vector.
                 do i = 1,n
                    vr(i) = eps3
                    vi(i) = zero
                 end do
              else
                 ! scale supplied initial vector.
                 norm = la_qlapy2(la_qnrm2(n,vr,1),la_qnrm2(n,vi,1))

                 rec = (eps3*rootn)/max(norm,nrmsml)
                 call la_qscal(n,rec,vr,1)
                 call la_qscal(n,rec,vi,1)
              end if
              if (rightv) then
                 ! lu decomposition with partial pivoting of b, replacing zero
                 ! pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(2,1) = -wi
                 do i = 2,n
                    b(i + 1,1) = zero
                 end do
                 loop_170: do i = 1,n - 1
                    absbii = la_qlapy2(b(i,i),b(i + 1,i))
                    ei = h(i + 1,i)
                    if (absbii < abs(ei)) then
                       ! interchange rows and eliminate.
                       xr = b(i,i)/ei
                       xi = b(i + 1,i)/ei
                       b(i,i) = ei
                       b(i + 1,i) = zero
                       do j = i + 1,n
                          temp = b(i + 1,j)
                          b(i + 1,j) = b(i,j) - xr*temp
                          b(j + 1,i + 1) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(i + 2,i) = -wi
                       b(i + 1,i + 1) = b(i + 1,i + 1) - xi*wi
                       b(i + 2,i + 1) = b(i + 2,i + 1) + xr*wi
                    else
                       ! eliminate without interchanging rows.
                       if (absbii == zero) then
                          b(i,i) = eps3
                          b(i + 1,i) = zero
                          absbii = eps3
                       end if
                       ei = (ei/absbii)/absbii
                       xr = b(i,i)*ei
                       xi = -b(i + 1,i)*ei
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j + 1,i + 1) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(i + 2,i + 1) = b(i + 2,i + 1) - wi
                    end if
                    ! compute 1-norm of offdiagonal elements of i-th row.
                    work(i) = la_qasum(n - i,b(i,i + 1),ldb) + la_qasum(n - i,b(i + 2, &
                              i),1)
                 end do loop_170
                 if (b(n,n) == zero .and. b(n + 1,n) == zero) b(n,n) = eps3
                 work(n) = zero
                 i1 = n
                 i2 = 1
                 i3 = -1
              else
                 ! ul decomposition with partial pivoting of conjg(b),
                 ! replacing zero pivots by eps3.
                 ! the imaginary part of the (i,j)-th element of u is stored in
                 ! b(j+1,i).
                 b(n + 1,n) = wi
                 do j = 1,n - 1
                    b(n + 1,j) = zero
                 end do
                 loop_210: do j = n,2,-1
                    ej = h(j,j - 1)
                    absbjj = la_qlapy2(b(j,j),b(j + 1,j))
                    if (absbjj < abs(ej)) then
                       ! interchange columns and eliminate
                       xr = b(j,j)/ej
                       xi = b(j + 1,j)/ej
                       b(j,j) = ej
                       b(j + 1,j) = zero
                       do i = 1,j - 1
                          temp = b(i,j - 1)
                          b(i,j - 1) = b(i,j) - xr*temp
                          b(j,i) = b(j + 1,i) - xi*temp
                          b(i,j) = temp
                          b(j + 1,i) = zero
                       end do
                       b(j + 1,j - 1) = wi
                       b(j - 1,j - 1) = b(j - 1,j - 1) + xi*wi
                       b(j,j - 1) = b(j,j - 1) - xr*wi
                    else
                       ! eliminate without interchange.
                       if (absbjj == zero) then
                          b(j,j) = eps3
                          b(j + 1,j) = zero
                          absbjj = eps3
                       end if
                       ej = (ej/absbjj)/absbjj
                       xr = b(j,j)*ej
                       xi = -b(j + 1,j)*ej
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - xr*b(i,j) + xi*b(j + 1,i)
                          b(j,i) = -xr*b(j + 1,i) - xi*b(i,j)
                       end do
                       b(j,j - 1) = b(j,j - 1) + wi
                    end if
                    ! compute 1-norm of offdiagonal elements of j-th column.
                    work(j) = la_qasum(j - 1,b(1,j),1) + la_qasum(j - 1,b(j + 1,1), &
                               ldb)
                 end do loop_210
                 if (b(1,1) == zero .and. b(2,1) == zero) b(1,1) = eps3
                 work(1) = zero
                 i1 = 1
                 i2 = n
                 i3 = 1
              end if
              loop_270: do its = 1,n
                 scale = one
                 vmax = one
                 vcrit = bignum
                 ! solve u*(xr,xi) = scale*(vr,vi) for a right eigenvector,
                   ! or u**t*(xr,xi) = scale*(vr,vi) for a left eigenvector,
                 ! overwriting (xr,xi) on (vr,vi).
                 loop_250: do i = i1,i2,i3
                    if (work(i) > vcrit) then
                       rec = one/vmax
                       call la_qscal(n,rec,vr,1)
                       call la_qscal(n,rec,vi,1)
                       scale = scale*rec
                       vmax = one
                       vcrit = bignum
                    end if
                    xr = vr(i)
                    xi = vi(i)
                    if (rightv) then
                       do j = i + 1,n
                          xr = xr - b(i,j)*vr(j) + b(j + 1,i)*vi(j)
                          xi = xi - b(i,j)*vi(j) - b(j + 1,i)*vr(j)
                       end do
                    else
                       do j = 1,i - 1
                          xr = xr - b(j,i)*vr(j) + b(i + 1,j)*vi(j)
                          xi = xi - b(j,i)*vi(j) - b(i + 1,j)*vr(j)
                       end do
                    end if
                    w = abs(b(i,i)) + abs(b(i + 1,i))
                    if (w > smlnum) then
                       if (w < one) then
                          w1 = abs(xr) + abs(xi)
                          if (w1 > w*bignum) then
                             rec = one/w1
                             call la_qscal(n,rec,vr,1)
                             call la_qscal(n,rec,vi,1)
                             xr = vr(i)
                             xi = vi(i)
                             scale = scale*rec
                             vmax = vmax*rec
                          end if
                       end if
                       ! divide by diagonal element of b.
                       call la_qladiv(xr,xi,b(i,i),b(i + 1,i),vr(i),vi(i))

                       vmax = max(abs(vr(i)) + abs(vi(i)),vmax)
                       vcrit = bignum/vmax
                    else
                       do j = 1,n
                          vr(j) = zero
                          vi(j) = zero
                       end do
                       vr(i) = one
                       vi(i) = one
                       scale = zero
                       vmax = one
                       vcrit = bignum
                    end if
                 end do loop_250
                 ! test for sufficient growth in the norm of (vr,vi).
                 vnorm = la_qasum(n,vr,1) + la_qasum(n,vi,1)
                 if (vnorm >= growto*scale) go to 280
                 ! choose a new orthogonal starting vector and try again.
                 y = eps3/(rootn + one)
                 vr(1) = eps3
                 vi(1) = zero
                 do i = 2,n
                    vr(i) = y
                    vi(i) = zero
                 end do
                 vr(n - its + 1) = vr(n - its + 1) - eps3*rootn
              end do loop_270
              ! failure to find eigenvector in n iterations
              info = 1
              280 continue
              ! normalize eigenvector.
              vnorm = zero
              do i = 1,n
                 vnorm = max(vnorm,abs(vr(i)) + abs(vi(i)))
              end do
              call la_qscal(n,one/vnorm,vr,1)
              call la_qscal(n,one/vnorm,vi,1)
           end if
           return
     end subroutine la_qlaein
#endif

     !> STREVC: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.

     pure subroutine la_strevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(sp),intent(in) :: t(ldt,*)
           real(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2
           real(sp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(sp) :: x(2,2)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_slamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_slabad(unfl,ovfl)
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
           n2 = 2*n
           if (rightv) then
              ! compute right eigenvectors.
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == 1) go to 130
                 if (ki == 1) go to 40
                 if (t(ki,ki - 1) == zero) go to 40
                 ip = -1
                 40 continue
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) go to 130
                    else
                       if (.not. select(ki - 1)) go to 130
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real right eigenvector
                    work(ki + n) = one
                    ! form right-hand side
                    do k = 1,ki - 1
                       work(k + n) = -t(k,ki)
                    end do
                    ! solve the upper quasi-triangular system:
                       ! (t(1:ki-1,1:ki-1) - wr)*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_slaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_sscal(ki,scale,work(1 + n),1)
                          work(j + n) = x(1,1)
                          ! update right-hand side
                          call la_saxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_slaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_sscal(ki,scale,work(1 + n),1)
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          ! update right-hand side
                          call la_saxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_saxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_scopy(ki,work(1 + n),1,vr(1,is),1)
                       ii = la_isamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_sscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 1) call la_sgemv('N',n,ki - 1,one,vr,ldvr,work(1 + n),1, &
                                 work(ki + n),vr(1,ki),1)
                       ii = la_isamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_sscal(n,remax,vr(1,ki),1)
                    end if
                 else
                    ! complex right eigenvector.
                    ! initial solve
                      ! [ (t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i* wi)]*x = 0.
                      ! [ (t(ki,ki-1)   t(ki,ki)   )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + n) = one
                       work(ki + n2) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + n) = -wi/t(ki,ki - 1)
                       work(ki + n2) = one
                    end if
                    work(ki + n) = zero
                    work(ki - 1 + n2) = zero
                    ! form right-hand side
                    do k = 1,ki - 2
                       work(k + n) = -work(ki - 1 + n)*t(k,ki - 1)
                       work(k + n2) = -work(ki + n2)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! (t(1:ki-2,1:ki-2) - (wr+i*wi))*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_slaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(ki,scale,work(1 + n),1)
                             call la_sscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          ! update the right-hand side
                          call la_saxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                          call la_saxpy(j - 1,-x(1,2),t(1,j),1,work(1 + n2),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_slaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(ki,scale,work(1 + n),1)
                             call la_sscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          work(j - 1 + n2) = x(1,2)
                          work(j + n2) = x(2,2)
                          ! update the right-hand side
                          call la_saxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_saxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                          call la_saxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + n2),1)

                          call la_saxpy(j - 2,-x(2,2),t(1,j),1,work(1 + n2),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_scopy(ki,work(1 + n),1,vr(1,is - 1),1)
                       call la_scopy(ki,work(1 + n2),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_sscal(ki,remax,vr(1,is - 1),1)
                       call la_sscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 2) then
                          call la_sgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n),1,work(ki - &
                                    1 + n),vr(1,ki - 1),1)
                          call la_sgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n2),1,work( &
                                    ki + n2),vr(1,ki),1)
                       else
                          call la_sscal(n,work(ki - 1 + n),vr(1,ki - 1),1)
                          call la_sscal(n,work(ki + n2),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_sscal(n,remax,vr(1,ki - 1),1)
                       call la_sscal(n,remax,vr(1,ki),1)
                    end if
                 end if
                 is = is - 1
                 if (ip /= 0) is = is - 1
                 130 continue
                 if (ip == 1) ip = 0
                 if (ip == -1) ip = 1
              end do loop_140
           end if
           if (leftv) then
              ! compute left eigenvectors.
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == -1) go to 250
                 if (ki == n) go to 150
                 if (t(ki + 1,ki) == zero) go to 150
                 ip = 1
                 150 continue
                 if (somev) then
                    if (.not. select(ki)) go to 250
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real left eigenvector.
                    work(ki + n) = one
                    ! form right-hand side
                    do k = ki + 1,n
                       work(k + n) = -t(ki,k)
                    end do
                    ! solve the quasi-triangular system:
                       ! (t(ki+1:n,ki+1:n) - wr)**t*x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_sdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          ! solve (t(j,j)-wr)**t*x = work
                          call la_slaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_sscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          vmax = max(abs(work(j + n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_sdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_sdot(j - ki - 1,t(ki + 1,j + 1),1, &
                                    work(ki + 1 + n),1)
                          ! solve
                            ! [t(j,j)-wr   t(j,j+1)     ]**t* x = scale*( work1 )
                            ! [t(j+1,j)    t(j+1,j+1)-wr]               ( work2 )
                          call la_slaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_sscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          work(j + 1 + n) = x(2,1)
                          vmax = max(abs(work(j + n)),abs(work(j + 1 + n)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_scopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       ii = la_isamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_sscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else
                       if (ki < n) call la_sgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + n),1,work(ki + n),vl(1,ki),1)
                       ii = la_isamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_sscal(n,remax,vl(1,ki),1)
                    end if
                 else
                    ! complex left eigenvector.
                     ! initial solve:
                       ! ((t(ki,ki)    t(ki,ki+1) )**t - (wr - i* wi))*x = 0.
                       ! ((t(ki+1,ki) t(ki+1,ki+1))                )
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + n) = wi/t(ki,ki + 1)
                       work(ki + 1 + n2) = one
                    else
                       work(ki + n) = one
                       work(ki + 1 + n2) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + n) = zero
                    work(ki + n2) = zero
                    ! form right-hand side
                    do k = ki + 2,n
                       work(k + n) = -work(ki + n)*t(ki,k)
                       work(k + n2) = -work(ki + 1 + n2)*t(ki + 1,k)
                    end do
                    ! solve complex quasi-triangular system:
                    ! ( t(ki+2,n:ki+2,n) - (wr-i*wi) )*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + n),1)
                             call la_sscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_sdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_sdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          ! solve (t(j,j)-(wr-i*wi))*(x11+i*x12)= wk+i*wk2
                          call la_slaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(n - ki + 1,scale,work(ki + n),1)
                             call la_sscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          vmax = max(abs(work(j + n)),abs(work(j + n2)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + n),1)
                             call la_sscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_sdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_sdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_sdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n),1)
                          work(j + 1 + n2) = work(j + 1 + n2) - la_sdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n2),1)
                          ! solve 2-by-2 complex linear equation
                            ! ([t(j,j)   t(j,j+1)  ]**t-(wr-i*wi)*i)*x = scale*b
                            ! ([t(j+1,j) t(j+1,j+1)]               )
                          call la_slaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(n - ki + 1,scale,work(ki + n),1)
                             call la_sscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          work(j + 1 + n) = x(2,1)
                          work(j + 1 + n2) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_scopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       call la_scopy(n - ki + 1,work(ki + n2),1,vl(ki,is + 1),1)
                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_sscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_sscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else
                       if (ki < n - 1) then
                          call la_sgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n),1,work(ki + n),vl(1,ki),1)
                          call la_sgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n2),1,work(ki + 1 + n2),vl(1,ki + 1),1)
                       else
                          call la_sscal(n,work(ki + n),vl(1,ki),1)
                          call la_sscal(n,work(ki + 1 + n2),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_sscal(n,remax,vl(1,ki),1)
                       call la_sscal(n,remax,vl(1,ki + 1),1)
                    end if
                 end if
                 is = is + 1
                 if (ip /= 0) is = is + 1
                 250 continue
                 if (ip == -1) ip = 0
                 if (ip == 1) ip = -1
              end do loop_260
           end if
           return
     end subroutine la_strevc
     !> DTREVC: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.

     pure subroutine la_dtrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(dp),intent(in) :: t(ldt,*)
           real(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2
           real(dp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(dp) :: x(2,2)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_dlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_dlabad(unfl,ovfl)
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
           n2 = 2*n
           if (rightv) then
              ! compute right eigenvectors.
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == 1) go to 130
                 if (ki == 1) go to 40
                 if (t(ki,ki - 1) == zero) go to 40
                 ip = -1
                 40 continue
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) go to 130
                    else
                       if (.not. select(ki - 1)) go to 130
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real right eigenvector
                    work(ki + n) = one
                    ! form right-hand side
                    do k = 1,ki - 1
                       work(k + n) = -t(k,ki)
                    end do
                    ! solve the upper quasi-triangular system:
                       ! (t(1:ki-1,1:ki-1) - wr)*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_dlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_dscal(ki,scale,work(1 + n),1)
                          work(j + n) = x(1,1)
                          ! update right-hand side
                          call la_daxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_dlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_dscal(ki,scale,work(1 + n),1)
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          ! update right-hand side
                          call la_daxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_daxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_dcopy(ki,work(1 + n),1,vr(1,is),1)
                       ii = la_idamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_dscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 1) call la_dgemv('N',n,ki - 1,one,vr,ldvr,work(1 + n),1, &
                                 work(ki + n),vr(1,ki),1)
                       ii = la_idamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_dscal(n,remax,vr(1,ki),1)
                    end if
                 else
                    ! complex right eigenvector.
                    ! initial solve
                      ! [ (t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i* wi)]*x = 0.
                      ! [ (t(ki,ki-1)   t(ki,ki)   )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + n) = one
                       work(ki + n2) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + n) = -wi/t(ki,ki - 1)
                       work(ki + n2) = one
                    end if
                    work(ki + n) = zero
                    work(ki - 1 + n2) = zero
                    ! form right-hand side
                    do k = 1,ki - 2
                       work(k + n) = -work(ki - 1 + n)*t(k,ki - 1)
                       work(k + n2) = -work(ki + n2)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! (t(1:ki-2,1:ki-2) - (wr+i*wi))*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_dlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(ki,scale,work(1 + n),1)
                             call la_dscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          ! update the right-hand side
                          call la_daxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                          call la_daxpy(j - 1,-x(1,2),t(1,j),1,work(1 + n2),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_dlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(ki,scale,work(1 + n),1)
                             call la_dscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          work(j - 1 + n2) = x(1,2)
                          work(j + n2) = x(2,2)
                          ! update the right-hand side
                          call la_daxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_daxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                          call la_daxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + n2),1)

                          call la_daxpy(j - 2,-x(2,2),t(1,j),1,work(1 + n2),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_dcopy(ki,work(1 + n),1,vr(1,is - 1),1)
                       call la_dcopy(ki,work(1 + n2),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_dscal(ki,remax,vr(1,is - 1),1)
                       call la_dscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 2) then
                          call la_dgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n),1,work(ki - &
                                    1 + n),vr(1,ki - 1),1)
                          call la_dgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n2),1,work( &
                                    ki + n2),vr(1,ki),1)
                       else
                          call la_dscal(n,work(ki - 1 + n),vr(1,ki - 1),1)
                          call la_dscal(n,work(ki + n2),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_dscal(n,remax,vr(1,ki - 1),1)
                       call la_dscal(n,remax,vr(1,ki),1)
                    end if
                 end if
                 is = is - 1
                 if (ip /= 0) is = is - 1
                 130 continue
                 if (ip == 1) ip = 0
                 if (ip == -1) ip = 1
              end do loop_140
           end if
           if (leftv) then
              ! compute left eigenvectors.
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == -1) go to 250
                 if (ki == n) go to 150
                 if (t(ki + 1,ki) == zero) go to 150
                 ip = 1
                 150 continue
                 if (somev) then
                    if (.not. select(ki)) go to 250
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real left eigenvector.
                    work(ki + n) = one
                    ! form right-hand side
                    do k = ki + 1,n
                       work(k + n) = -t(ki,k)
                    end do
                    ! solve the quasi-triangular system:
                       ! (t(ki+1:n,ki+1:n) - wr)**t*x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_ddot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          ! solve (t(j,j)-wr)**t*x = work
                          call la_dlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_dscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          vmax = max(abs(work(j + n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_ddot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_ddot(j - ki - 1,t(ki + 1,j + 1),1, &
                                    work(ki + 1 + n),1)
                          ! solve
                            ! [t(j,j)-wr   t(j,j+1)     ]**t * x = scale*( work1 )
                            ! [t(j+1,j)    t(j+1,j+1)-wr]                ( work2 )
                          call la_dlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_dscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          work(j + 1 + n) = x(2,1)
                          vmax = max(abs(work(j + n)),abs(work(j + 1 + n)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_dcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       ii = la_idamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_dscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else
                       if (ki < n) call la_dgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + n),1,work(ki + n),vl(1,ki),1)
                       ii = la_idamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_dscal(n,remax,vl(1,ki),1)
                    end if
                 else
                    ! complex left eigenvector.
                     ! initial solve:
                       ! ((t(ki,ki)    t(ki,ki+1) )**t - (wr - i* wi))*x = 0.
                       ! ((t(ki+1,ki) t(ki+1,ki+1))                )
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + n) = wi/t(ki,ki + 1)
                       work(ki + 1 + n2) = one
                    else
                       work(ki + n) = one
                       work(ki + 1 + n2) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + n) = zero
                    work(ki + n2) = zero
                    ! form right-hand side
                    do k = ki + 2,n
                       work(k + n) = -work(ki + n)*t(ki,k)
                       work(k + n2) = -work(ki + 1 + n2)*t(ki + 1,k)
                    end do
                    ! solve complex quasi-triangular system:
                    ! ( t(ki+2,n:ki+2,n) - (wr-i*wi) )*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + n),1)
                             call la_dscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_ddot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_ddot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          ! solve (t(j,j)-(wr-i*wi))*(x11+i*x12)= wk+i*wk2
                          call la_dlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(n - ki + 1,scale,work(ki + n),1)
                             call la_dscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          vmax = max(abs(work(j + n)),abs(work(j + n2)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + n),1)
                             call la_dscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_ddot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_ddot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_ddot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n),1)
                          work(j + 1 + n2) = work(j + 1 + n2) - la_ddot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n2),1)
                          ! solve 2-by-2 complex linear equation
                            ! ([t(j,j)   t(j,j+1)  ]**t-(wr-i*wi)*i)*x = scale*b
                            ! ([t(j+1,j) t(j+1,j+1)]               )
                          call la_dlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(n - ki + 1,scale,work(ki + n),1)
                             call la_dscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          work(j + 1 + n) = x(2,1)
                          work(j + 1 + n2) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_dcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       call la_dcopy(n - ki + 1,work(ki + n2),1,vl(ki,is + 1),1)
                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_dscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_dscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else
                       if (ki < n - 1) then
                          call la_dgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n),1,work(ki + n),vl(1,ki),1)
                          call la_dgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n2),1,work(ki + 1 + n2),vl(1,ki + 1),1)
                       else
                          call la_dscal(n,work(ki + n),vl(1,ki),1)
                          call la_dscal(n,work(ki + 1 + n2),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_dscal(n,remax,vl(1,ki),1)
                       call la_dscal(n,remax,vl(1,ki + 1),1)
                    end if
                 end if
                 is = is + 1
                 if (ip /= 0) is = is + 1
                 250 continue
                 if (ip == -1) ip = 0
                 if (ip == 1) ip = -1
              end do loop_260
           end if
           return
     end subroutine la_dtrevc
#ifdef LA_WITH_XDP
     !> XTREVC: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by XHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.

     pure subroutine la_xtrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(xdp),intent(in) :: t(ldt,*)
           real(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2
           real(xdp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(xdp) :: x(2,2)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_xlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_xlabad(unfl,ovfl)
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
           n2 = 2*n
           if (rightv) then
              ! compute right eigenvectors.
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == 1) go to 130
                 if (ki == 1) go to 40
                 if (t(ki,ki - 1) == zero) go to 40
                 ip = -1
                 40 continue
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) go to 130
                    else
                       if (.not. select(ki - 1)) go to 130
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real right eigenvector
                    work(ki + n) = one
                    ! form right-hand side
                    do k = 1,ki - 1
                       work(k + n) = -t(k,ki)
                    end do
                    ! solve the upper quasi-triangular system:
                       ! (t(1:ki-1,1:ki-1) - wr)*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_xlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_xscal(ki,scale,work(1 + n),1)
                          work(j + n) = x(1,1)
                          ! update right-hand side
                          call la_xaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_xlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_xscal(ki,scale,work(1 + n),1)
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          ! update right-hand side
                          call la_xaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_xaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_xcopy(ki,work(1 + n),1,vr(1,is),1)
                       ii = la_ixamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_xscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 1) call la_xgemv('N',n,ki - 1,one,vr,ldvr,work(1 + n),1, &
                                 work(ki + n),vr(1,ki),1)
                       ii = la_ixamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_xscal(n,remax,vr(1,ki),1)
                    end if
                 else
                    ! complex right eigenvector.
                    ! initial solve
                      ! [ (t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i* wi)]*x = 0.
                      ! [ (t(ki,ki-1)   t(ki,ki)   )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + n) = one
                       work(ki + n2) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + n) = -wi/t(ki,ki - 1)
                       work(ki + n2) = one
                    end if
                    work(ki + n) = zero
                    work(ki - 1 + n2) = zero
                    ! form right-hand side
                    do k = 1,ki - 2
                       work(k + n) = -work(ki - 1 + n)*t(k,ki - 1)
                       work(k + n2) = -work(ki + n2)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! (t(1:ki-2,1:ki-2) - (wr+i*wi))*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_xlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(ki,scale,work(1 + n),1)
                             call la_xscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          ! update the right-hand side
                          call la_xaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                          call la_xaxpy(j - 1,-x(1,2),t(1,j),1,work(1 + n2),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_xlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(ki,scale,work(1 + n),1)
                             call la_xscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          work(j - 1 + n2) = x(1,2)
                          work(j + n2) = x(2,2)
                          ! update the right-hand side
                          call la_xaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_xaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                          call la_xaxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + n2),1)

                          call la_xaxpy(j - 2,-x(2,2),t(1,j),1,work(1 + n2),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_xcopy(ki,work(1 + n),1,vr(1,is - 1),1)
                       call la_xcopy(ki,work(1 + n2),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_xscal(ki,remax,vr(1,is - 1),1)
                       call la_xscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 2) then
                          call la_xgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n),1,work(ki - &
                                    1 + n),vr(1,ki - 1),1)
                          call la_xgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n2),1,work( &
                                    ki + n2),vr(1,ki),1)
                       else
                          call la_xscal(n,work(ki - 1 + n),vr(1,ki - 1),1)
                          call la_xscal(n,work(ki + n2),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_xscal(n,remax,vr(1,ki - 1),1)
                       call la_xscal(n,remax,vr(1,ki),1)
                    end if
                 end if
                 is = is - 1
                 if (ip /= 0) is = is - 1
                 130 continue
                 if (ip == 1) ip = 0
                 if (ip == -1) ip = 1
              end do loop_140
           end if
           if (leftv) then
              ! compute left eigenvectors.
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == -1) go to 250
                 if (ki == n) go to 150
                 if (t(ki + 1,ki) == zero) go to 150
                 ip = 1
                 150 continue
                 if (somev) then
                    if (.not. select(ki)) go to 250
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real left eigenvector.
                    work(ki + n) = one
                    ! form right-hand side
                    do k = ki + 1,n
                       work(k + n) = -t(ki,k)
                    end do
                    ! solve the quasi-triangular system:
                       ! (t(ki+1:n,ki+1:n) - wr)**t*x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_xdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          ! solve (t(j,j)-wr)**t*x = work
                          call la_xlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_xscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          vmax = max(abs(work(j + n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_xdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_xdot(j - ki - 1,t(ki + 1,j + 1),1, &
                                    work(ki + 1 + n),1)
                          ! solve
                            ! [t(j,j)-wr   t(j,j+1)     ]**t * x = scale*( work1 )
                            ! [t(j+1,j)    t(j+1,j+1)-wr]                ( work2 )
                          call la_xlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_xscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          work(j + 1 + n) = x(2,1)
                          vmax = max(abs(work(j + n)),abs(work(j + 1 + n)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_xcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       ii = la_ixamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_xscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else
                       if (ki < n) call la_xgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + n),1,work(ki + n),vl(1,ki),1)
                       ii = la_ixamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_xscal(n,remax,vl(1,ki),1)
                    end if
                 else
                    ! complex left eigenvector.
                     ! initial solve:
                       ! ((t(ki,ki)    t(ki,ki+1) )**t - (wr - i* wi))*x = 0.
                       ! ((t(ki+1,ki) t(ki+1,ki+1))                )
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + n) = wi/t(ki,ki + 1)
                       work(ki + 1 + n2) = one
                    else
                       work(ki + n) = one
                       work(ki + 1 + n2) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + n) = zero
                    work(ki + n2) = zero
                    ! form right-hand side
                    do k = ki + 2,n
                       work(k + n) = -work(ki + n)*t(ki,k)
                       work(k + n2) = -work(ki + 1 + n2)*t(ki + 1,k)
                    end do
                    ! solve complex quasi-triangular system:
                    ! ( t(ki+2,n:ki+2,n) - (wr-i*wi) )*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + n),1)
                             call la_xscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_xdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_xdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          ! solve (t(j,j)-(wr-i*wi))*(x11+i*x12)= wk+i*wk2
                          call la_xlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(n - ki + 1,scale,work(ki + n),1)
                             call la_xscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          vmax = max(abs(work(j + n)),abs(work(j + n2)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + n),1)
                             call la_xscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_xdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_xdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_xdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n),1)
                          work(j + 1 + n2) = work(j + 1 + n2) - la_xdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n2),1)
                          ! solve 2-by-2 complex linear equation
                            ! ([t(j,j)   t(j,j+1)  ]**t-(wr-i*wi)*i)*x = scale*b
                            ! ([t(j+1,j) t(j+1,j+1)]               )
                          call la_xlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(n - ki + 1,scale,work(ki + n),1)
                             call la_xscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          work(j + 1 + n) = x(2,1)
                          work(j + 1 + n2) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_xcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       call la_xcopy(n - ki + 1,work(ki + n2),1,vl(ki,is + 1),1)
                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_xscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_xscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else
                       if (ki < n - 1) then
                          call la_xgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n),1,work(ki + n),vl(1,ki),1)
                          call la_xgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n2),1,work(ki + 1 + n2),vl(1,ki + 1),1)
                       else
                          call la_xscal(n,work(ki + n),vl(1,ki),1)
                          call la_xscal(n,work(ki + 1 + n2),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_xscal(n,remax,vl(1,ki),1)
                       call la_xscal(n,remax,vl(1,ki + 1),1)
                    end if
                 end if
                 is = is + 1
                 if (ip /= 0) is = is + 1
                 250 continue
                 if (ip == -1) ip = 0
                 if (ip == 1) ip = -1
              end do loop_260
           end if
           return
     end subroutine la_xtrevc
#endif
#ifdef LA_WITH_QP
     !> QTREVC: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by QHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.

     pure subroutine la_qtrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(qp),intent(in) :: t(ldt,*)
           real(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,n2
           real(qp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(qp) :: x(2,2)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_qlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_qlabad(unfl,ovfl)
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
           n2 = 2*n
           if (rightv) then
              ! compute right eigenvectors.
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == 1) go to 130
                 if (ki == 1) go to 40
                 if (t(ki,ki - 1) == zero) go to 40
                 ip = -1
                 40 continue
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) go to 130
                    else
                       if (.not. select(ki - 1)) go to 130
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real right eigenvector
                    work(ki + n) = one
                    ! form right-hand side
                    do k = 1,ki - 1
                       work(k + n) = -t(k,ki)
                    end do
                    ! solve the upper quasi-triangular system:
                       ! (t(1:ki-1,1:ki-1) - wr)*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_qlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_qscal(ki,scale,work(1 + n),1)
                          work(j + n) = x(1,1)
                          ! update right-hand side
                          call la_qaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_qlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_qscal(ki,scale,work(1 + n),1)
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          ! update right-hand side
                          call la_qaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_qaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_qcopy(ki,work(1 + n),1,vr(1,is),1)
                       ii = la_iqamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_qscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 1) call la_qgemv('N',n,ki - 1,one,vr,ldvr,work(1 + n),1, &
                                 work(ki + n),vr(1,ki),1)
                       ii = la_iqamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_qscal(n,remax,vr(1,ki),1)
                    end if
                 else
                    ! complex right eigenvector.
                    ! initial solve
                      ! [ (t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i* wi)]*x = 0.
                      ! [ (t(ki,ki-1)   t(ki,ki)   )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + n) = one
                       work(ki + n2) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + n) = -wi/t(ki,ki - 1)
                       work(ki + n2) = one
                    end if
                    work(ki + n) = zero
                    work(ki - 1 + n2) = zero
                    ! form right-hand side
                    do k = 1,ki - 2
                       work(k + n) = -work(ki - 1 + n)*t(k,ki - 1)
                       work(k + n2) = -work(ki + n2)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! (t(1:ki-2,1:ki-2) - (wr+i*wi))*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_qlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(ki,scale,work(1 + n),1)
                             call la_qscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          ! update the right-hand side
                          call la_qaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + n),1)

                          call la_qaxpy(j - 1,-x(1,2),t(1,j),1,work(1 + n2),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_qlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(ki,scale,work(1 + n),1)
                             call la_qscal(ki,scale,work(1 + n2),1)
                          end if
                          work(j - 1 + n) = x(1,1)
                          work(j + n) = x(2,1)
                          work(j - 1 + n2) = x(1,2)
                          work(j + n2) = x(2,2)
                          ! update the right-hand side
                          call la_qaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + n),1)

                          call la_qaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + n),1)

                          call la_qaxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + n2),1)

                          call la_qaxpy(j - 2,-x(2,2),t(1,j),1,work(1 + n2),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       call la_qcopy(ki,work(1 + n),1,vr(1,is - 1),1)
                       call la_qcopy(ki,work(1 + n2),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_qscal(ki,remax,vr(1,is - 1),1)
                       call la_qscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else
                       if (ki > 2) then
                          call la_qgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n),1,work(ki - &
                                    1 + n),vr(1,ki - 1),1)
                          call la_qgemv('N',n,ki - 2,one,vr,ldvr,work(1 + n2),1,work( &
                                    ki + n2),vr(1,ki),1)
                       else
                          call la_qscal(n,work(ki - 1 + n),vr(1,ki - 1),1)
                          call la_qscal(n,work(ki + n2),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_qscal(n,remax,vr(1,ki - 1),1)
                       call la_qscal(n,remax,vr(1,ki),1)
                    end if
                 end if
                 is = is - 1
                 if (ip /= 0) is = is - 1
                 130 continue
                 if (ip == 1) ip = 0
                 if (ip == -1) ip = 1
              end do loop_140
           end if
           if (leftv) then
              ! compute left eigenvectors.
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == -1) go to 250
                 if (ki == n) go to 150
                 if (t(ki + 1,ki) == zero) go to 150
                 ip = 1
                 150 continue
                 if (somev) then
                    if (.not. select(ki)) go to 250
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! real left eigenvector.
                    work(ki + n) = one
                    ! form right-hand side
                    do k = ki + 1,n
                       work(k + n) = -t(ki,k)
                    end do
                    ! solve the quasi-triangular system:
                       ! (t(ki+1:n,ki+1:n) - wr)**t*x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_qdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          ! solve (t(j,j)-wr)**t*x = work
                          call la_qlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_qscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          vmax = max(abs(work(j + n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_qdot(j - ki - 1,t(ki + 1,j),1,work( &
                                    ki + 1 + n),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_qdot(j - ki - 1,t(ki + 1,j + 1),1, &
                                    work(ki + 1 + n),1)
                          ! solve
                            ! [t(j,j)-wr   t(j,j+1)     ]**t * x = scale*( work1 )
                            ! [t(j+1,j)    t(j+1,j+1)-wr]                ( work2 )
                          call la_qlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_qscal(n - ki + 1,scale,work(ki + n),1)

                          work(j + n) = x(1,1)
                          work(j + 1 + n) = x(2,1)
                          vmax = max(abs(work(j + n)),abs(work(j + 1 + n)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_qcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       ii = la_iqamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_qscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else
                       if (ki < n) call la_qgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + n),1,work(ki + n),vl(1,ki),1)
                       ii = la_iqamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_qscal(n,remax,vl(1,ki),1)
                    end if
                 else
                    ! complex left eigenvector.
                     ! initial solve:
                       ! ((t(ki,ki)    t(ki,ki+1) )**t - (wr - i* wi))*x = 0.
                       ! ((t(ki+1,ki) t(ki+1,ki+1))                )
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + n) = wi/t(ki,ki + 1)
                       work(ki + 1 + n2) = one
                    else
                       work(ki + n) = one
                       work(ki + 1 + n2) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + n) = zero
                    work(ki + n2) = zero
                    ! form right-hand side
                    do k = ki + 2,n
                       work(k + n) = -work(ki + n)*t(ki,k)
                       work(k + n2) = -work(ki + 1 + n2)*t(ki + 1,k)
                    end do
                    ! solve complex quasi-triangular system:
                    ! ( t(ki+2,n:ki+2,n) - (wr-i*wi) )*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + n),1)
                             call la_qscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_qdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_qdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          ! solve (t(j,j)-(wr-i*wi))*(x11+i*x12)= wk+i*wk2
                          call la_qlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(n - ki + 1,scale,work(ki + n),1)
                             call la_qscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          vmax = max(abs(work(j + n)),abs(work(j + n2)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + n),1)
                             call la_qscal(n - ki + 1,rec,work(ki + n2),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + n) = work(j + n) - la_qdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n),1)
                          work(j + n2) = work(j + n2) - la_qdot(j - ki - 2,t(ki + 2,j),1,work( &
                                    ki + 2 + n2),1)
                          work(j + 1 + n) = work(j + 1 + n) - la_qdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n),1)
                          work(j + 1 + n2) = work(j + 1 + n2) - la_qdot(j - ki - 2,t(ki + 2,j + 1),1, &
                                    work(ki + 2 + n2),1)
                          ! solve 2-by-2 complex linear equation
                            ! ([t(j,j)   t(j,j+1)  ]**t-(wr-i*wi)*i)*x = scale*b
                            ! ([t(j+1,j) t(j+1,j+1)]               )
                          call la_qlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(n - ki + 1,scale,work(ki + n),1)
                             call la_qscal(n - ki + 1,scale,work(ki + n2),1)
                          end if
                          work(j + n) = x(1,1)
                          work(j + n2) = x(1,2)
                          work(j + 1 + n) = x(2,1)
                          work(j + 1 + n2) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       call la_qcopy(n - ki + 1,work(ki + n),1,vl(ki,is),1)
                       call la_qcopy(n - ki + 1,work(ki + n2),1,vl(ki,is + 1),1)
                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_qscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_qscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else
                       if (ki < n - 1) then
                          call la_qgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n),1,work(ki + n),vl(1,ki),1)
                          call la_qgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    n2),1,work(ki + 1 + n2),vl(1,ki + 1),1)
                       else
                          call la_qscal(n,work(ki + n),vl(1,ki),1)
                          call la_qscal(n,work(ki + 1 + n2),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_qscal(n,remax,vl(1,ki),1)
                       call la_qscal(n,remax,vl(1,ki + 1),1)
                    end if
                 end if
                 is = is + 1
                 if (ip /= 0) is = is + 1
                 250 continue
                 if (ip == -1) ip = 0
                 if (ip == 1) ip = -1
              end do loop_260
           end if
           return
     end subroutine la_qtrevc
#endif

     !> STREVC3: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**T)*T = w*(y**T)
     !> where y**T denotes the transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_strevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(sp),intent(in) :: t(ldt,*)
           real(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,iv,maxwrk,nb, &
                     ki2
           real(sp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(sp) :: x(2,2)
           integer(ilp) :: iscomplex(nbmax)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           nb = la_ilaenv(1,'STREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           lquery = (lwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -14
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_slaset('F',n,1 + 2*nb,zero,zero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_slamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_slabad(unfl,ovfl)
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first  of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
             ! iscomplex array stores ip for each column in current block.
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! for complex right vector, uses iv-1 for real part and iv for complex part.
              ! non-blocked version always uses iv=2;
              ! blocked     version starts with iv=nb, goes down to 1 or 2.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 2
              if (nb > 2) then
                 iv = nb
              end if
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == -1) then
                    ! previous iteration (ki+1) was second of conjugate pair,
                    ! so this ki is first of conjugate pair; skip to end of loop
                    ip = 1
                    cycle loop_140
                 else if (ki == 1) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki,ki - 1) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is second of conjugate pair
                    ip = -1
                 end if
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) cycle loop_140
                    else
                       if (.not. select(ki - 1)) cycle loop_140
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real right eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = 1,ki - 1
                       work(k + iv*n) = -t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-1,1:ki-1) - wr ]*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_slaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_sscal(ki,scale,work(1 + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          ! update right-hand side
                          call la_saxpy(j - 1,-x(1,1),t(1,j),1,work(1 + iv*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_slaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_sscal(ki,scale,work(1 + iv*n),1)

                          work(j - 1 + iv*n) = x(1,1)
                          work(j + iv*n) = x(2,1)
                          ! update right-hand side
                          call la_saxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + iv*n),1)

                          call la_saxpy(j - 2,-x(2,1),t(1,j),1,work(1 + iv*n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_scopy(ki,work(1 + iv*n),1,vr(1,is),1)
                       ii = la_isamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_sscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 1) call la_sgemv('N',n,ki - 1,one,vr,ldvr,work(1 + iv*n), &
                                 1,work(ki + iv*n),vr(1,ki),1)
                       ii = la_isamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_sscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex right eigenvector.
                    ! initial solve
                    ! [ ( t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i*wi) ]*x = 0.
                    ! [ ( t(ki,  ki-1) t(ki,  ki) )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + (iv - 1)*n) = one
                       work(ki + (iv)*n) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + (iv - 1)*n) = -wi/t(ki,ki - 1)
                       work(ki + (iv)*n) = one
                    end if
                    work(ki + (iv - 1)*n) = zero
                    work(ki - 1 + (iv)*n) = zero
                    ! form right-hand side.
                    do k = 1,ki - 2
                       work(k + (iv - 1)*n) = -work(ki - 1 + (iv - 1)*n)*t(k,ki - 1)
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-2,1:ki-2) - (wr+i*wi) ]*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_slaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_sscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j + (iv - 1)*n) = x(1,1)
                          work(j + (iv)*n) = x(1,2)
                          ! update the right-hand side
                          call la_saxpy(j - 1,-x(1,1),t(1,j),1,work(1 + (iv - 1)*n),1)

                          call la_saxpy(j - 1,-x(1,2),t(1,j),1,work(1 + (iv)*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_slaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_sscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j - 1 + (iv - 1)*n) = x(1,1)
                          work(j + (iv - 1)*n) = x(2,1)
                          work(j - 1 + (iv)*n) = x(1,2)
                          work(j + (iv)*n) = x(2,2)
                          ! update the right-hand side
                          call la_saxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + (iv - 1)*n), &
                                     1)
                          call la_saxpy(j - 2,-x(2,1),t(1,j),1,work(1 + (iv - 1)*n), &
                                    1)
                          call la_saxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + (iv)*n), &
                                    1)
                          call la_saxpy(j - 2,-x(2,2),t(1,j),1,work(1 + (iv)*n),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_scopy(ki,work(1 + (iv - 1)*n),1,vr(1,is - 1),1)
                       call la_scopy(ki,work(1 + (iv)*n),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_sscal(ki,remax,vr(1,is - 1),1)
                       call la_sscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 2) then
                          call la_sgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv - 1)*n), &
                                    1,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_sgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv)*n),1, &
                                    work(ki + (iv)*n),vr(1,ki),1)
                       else
                          call la_sscal(n,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_sscal(n,work(ki + (iv)*n),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_sscal(n,remax,vr(1,ki - 1),1)
                       call la_sscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + (iv - 1)*n) = zero
                          work(k + (iv)*n) = zero
                       end do
                       iscomplex(iv - 1) = -ip
                       iscomplex(iv) = ip
                       iv = iv - 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki-1 and ki)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki - 1
                    end if
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv <= 2) .or. (ki2 == 1)) then
                       call la_sgemm('N','N',n,nb - iv + 1,ki2 + nb - iv,one,vr,ldvr,work(1 + &
                                 (iv)*n),n,zero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_isamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_sscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_slacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki2), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if ! blocked back-transform
                 is = is - 1
                 if (ip /= 0) is = is - 1
              end do loop_140
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! for complex left vector, uses iv for real part and iv+1 for complex part.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb-1 or nb.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 1
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == 1) then
                    ! previous iteration (ki-1) was first of conjugate pair,
                    ! so this ki is second of conjugate pair; skip to end of loop
                    ip = -1
                    cycle loop_260
                 else if (ki == n) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki + 1,ki) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is first of conjugate pair
                    ip = 1
                 end if
                 if (somev) then
                    if (.not. select(ki)) cycle loop_260
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real left eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = ki + 1,n
                       work(k + iv*n) = -t(ki,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+1:n,ki+1:n) - wr ]**t * x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_sdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          ! solve [ t(j,j) - wr ]**t * x = work
                          call la_slaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_sscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          vmax = max(abs(work(j + iv*n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_sdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          work(j + 1 + iv*n) = work(j + 1 + iv*n) - la_sdot(j - ki - 1,t(ki + 1,j + 1) &
                                    ,1,work(ki + 1 + iv*n),1)
                          ! solve
                          ! [ t(j,j)-wr   t(j,j+1)      ]**t * x = scale*( work1 )
                          ! [ t(j+1,j)    t(j+1,j+1)-wr ]                ( work2 )
                          call la_slaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_sscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          work(j + 1 + iv*n) = x(2,1)
                          vmax = max(abs(work(j + iv*n)),abs(work(j + 1 + iv*n)),vmax)

                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_scopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                       ii = la_isamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_sscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n) call la_sgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + iv*n),1,work(ki + iv*n),vl(1,ki),1)
                       ii = la_isamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_sscal(n,remax,vl(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex left eigenvector.
                    ! initial solve:
                    ! [ ( t(ki,ki)    t(ki,ki+1)  )**t - (wr - i* wi) ]*x = 0.
                    ! [ ( t(ki+1,ki) t(ki+1,ki+1) )                   ]
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + (iv)*n) = wi/t(ki,ki + 1)
                       work(ki + 1 + (iv + 1)*n) = one
                    else
                       work(ki + (iv)*n) = one
                       work(ki + 1 + (iv + 1)*n) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + (iv)*n) = zero
                    work(ki + (iv + 1)*n) = zero
                    ! form right-hand side.
                    do k = ki + 2,n
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(ki,k)
                       work(k + (iv + 1)*n) = -work(ki + 1 + (iv + 1)*n)*t(ki + 1,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+2:n,ki+2:n)**t - (wr-i*wi) ]*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_sscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_sdot(j - ki - 2,t(ki + 2,j) &
                                    ,1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_sdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve [ t(j,j)-(wr-i*wi) ]*(x11+i*x12)= wk+i*wk2
                          call la_slaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_sscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          vmax = max(abs(work(j + (iv)*n)),abs(work(j + (iv + 1)*n)),vmax)

                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_sscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_sscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_sdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_sdot(j - ki - 2,t(ki + 2, &
                                     j),1,work(ki + 2 + (iv + 1)*n),1)
                          work(j + 1 + (iv)*n) = work(j + 1 + (iv)*n) - la_sdot(j - ki - 2,t(ki + 2, &
                                     j + 1),1,work(ki + 2 + (iv)*n),1)
                          work(j + 1 + (iv + 1)*n) = work(j + 1 + (iv + 1)*n) - la_sdot(j - ki - 2,t(ki + &
                                    2,j + 1),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve 2-by-2 complex linear equation
                          ! [ (t(j,j)   t(j,j+1)  )**t - (wr-i*wi)*i ]*x = scale*b
                          ! [ (t(j+1,j) t(j+1,j+1))                  ]
                          call la_slaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_sscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_sscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          work(j + 1 + (iv)*n) = x(2,1)
                          work(j + 1 + (iv + 1)*n) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_scopy(n - ki + 1,work(ki + (iv)*n),1,vl(ki,is),1)

                       call la_scopy(n - ki + 1,work(ki + (iv + 1)*n),1,vl(ki,is + 1),1)

                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_sscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_sscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n - 1) then
                          call la_sgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv)*n),1,work(ki + (iv)*n),vl(1,ki),1)
                          call la_sgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv + 1)*n),1,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       else
                          call la_sscal(n,work(ki + (iv)*n),vl(1,ki),1)
                          call la_sscal(n,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_sscal(n,remax,vl(1,ki),1)
                       call la_sscal(n,remax,vl(1,ki + 1),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + (iv)*n) = zero
                          work(k + (iv + 1)*n) = zero
                       end do
                       iscomplex(iv) = ip
                       iscomplex(iv + 1) = -ip
                       iv = iv + 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki and ki+1)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki + 1
                    end if
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv >= nb - 1) .or. (ki2 == n)) then
                       call la_sgemm('N','N',n,iv,n - ki2 + iv,one,vl(1,ki2 - iv + 1),ldvl, &
                                 work(ki2 - iv + 1 + (1)*n),n,zero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_isamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_sscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_slacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki2 - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if ! blocked back-transform
                 is = is + 1
                 if (ip /= 0) is = is + 1
              end do loop_260
           end if
           return
     end subroutine la_strevc3
     !> DTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**T)*T = w*(y**T)
     !> where y**T denotes the transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_dtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(dp),intent(in) :: t(ldt,*)
           real(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,iv,maxwrk,nb, &
                     ki2
           real(dp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(dp) :: x(2,2)
           integer(ilp) :: iscomplex(nbmax)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           nb = la_ilaenv(1,'DTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           lquery = (lwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -14
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_dlaset('F',n,1 + 2*nb,zero,zero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_dlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_dlabad(unfl,ovfl)
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first  of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
             ! iscomplex array stores ip for each column in current block.
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! for complex right vector, uses iv-1 for real part and iv for complex part.
              ! non-blocked version always uses iv=2;
              ! blocked     version starts with iv=nb, goes down to 1 or 2.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 2
              if (nb > 2) then
                 iv = nb
              end if
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == -1) then
                    ! previous iteration (ki+1) was second of conjugate pair,
                    ! so this ki is first of conjugate pair; skip to end of loop
                    ip = 1
                    cycle loop_140
                 else if (ki == 1) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki,ki - 1) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is second of conjugate pair
                    ip = -1
                 end if
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) cycle loop_140
                    else
                       if (.not. select(ki - 1)) cycle loop_140
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real right eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = 1,ki - 1
                       work(k + iv*n) = -t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-1,1:ki-1) - wr ]*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_dlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_dscal(ki,scale,work(1 + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          ! update right-hand side
                          call la_daxpy(j - 1,-x(1,1),t(1,j),1,work(1 + iv*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_dlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_dscal(ki,scale,work(1 + iv*n),1)

                          work(j - 1 + iv*n) = x(1,1)
                          work(j + iv*n) = x(2,1)
                          ! update right-hand side
                          call la_daxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + iv*n),1)

                          call la_daxpy(j - 2,-x(2,1),t(1,j),1,work(1 + iv*n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_dcopy(ki,work(1 + iv*n),1,vr(1,is),1)
                       ii = la_idamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_dscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 1) call la_dgemv('N',n,ki - 1,one,vr,ldvr,work(1 + iv*n), &
                                 1,work(ki + iv*n),vr(1,ki),1)
                       ii = la_idamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_dscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex right eigenvector.
                    ! initial solve
                    ! [ ( t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i*wi) ]*x = 0.
                    ! [ ( t(ki,  ki-1) t(ki,  ki) )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + (iv - 1)*n) = one
                       work(ki + (iv)*n) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + (iv - 1)*n) = -wi/t(ki,ki - 1)
                       work(ki + (iv)*n) = one
                    end if
                    work(ki + (iv - 1)*n) = zero
                    work(ki - 1 + (iv)*n) = zero
                    ! form right-hand side.
                    do k = 1,ki - 2
                       work(k + (iv - 1)*n) = -work(ki - 1 + (iv - 1)*n)*t(k,ki - 1)
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-2,1:ki-2) - (wr+i*wi) ]*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_dlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_dscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j + (iv - 1)*n) = x(1,1)
                          work(j + (iv)*n) = x(1,2)
                          ! update the right-hand side
                          call la_daxpy(j - 1,-x(1,1),t(1,j),1,work(1 + (iv - 1)*n),1)

                          call la_daxpy(j - 1,-x(1,2),t(1,j),1,work(1 + (iv)*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_dlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_dscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j - 1 + (iv - 1)*n) = x(1,1)
                          work(j + (iv - 1)*n) = x(2,1)
                          work(j - 1 + (iv)*n) = x(1,2)
                          work(j + (iv)*n) = x(2,2)
                          ! update the right-hand side
                          call la_daxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + (iv - 1)*n), &
                                     1)
                          call la_daxpy(j - 2,-x(2,1),t(1,j),1,work(1 + (iv - 1)*n), &
                                    1)
                          call la_daxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + (iv)*n), &
                                    1)
                          call la_daxpy(j - 2,-x(2,2),t(1,j),1,work(1 + (iv)*n),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_dcopy(ki,work(1 + (iv - 1)*n),1,vr(1,is - 1),1)
                       call la_dcopy(ki,work(1 + (iv)*n),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_dscal(ki,remax,vr(1,is - 1),1)
                       call la_dscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 2) then
                          call la_dgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv - 1)*n), &
                                    1,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_dgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv)*n),1, &
                                    work(ki + (iv)*n),vr(1,ki),1)
                       else
                          call la_dscal(n,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_dscal(n,work(ki + (iv)*n),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_dscal(n,remax,vr(1,ki - 1),1)
                       call la_dscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + (iv - 1)*n) = zero
                          work(k + (iv)*n) = zero
                       end do
                       iscomplex(iv - 1) = -ip
                       iscomplex(iv) = ip
                       iv = iv - 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki-1 and ki)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki - 1
                    end if
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv <= 2) .or. (ki2 == 1)) then
                       call la_dgemm('N','N',n,nb - iv + 1,ki2 + nb - iv,one,vr,ldvr,work(1 + &
                                 (iv)*n),n,zero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_idamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_dscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_dlacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki2), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if ! blocked back-transform
                 is = is - 1
                 if (ip /= 0) is = is - 1
              end do loop_140
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! for complex left vector, uses iv for real part and iv+1 for complex part.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb-1 or nb.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 1
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == 1) then
                    ! previous iteration (ki-1) was first of conjugate pair,
                    ! so this ki is second of conjugate pair; skip to end of loop
                    ip = -1
                    cycle loop_260
                 else if (ki == n) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki + 1,ki) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is first of conjugate pair
                    ip = 1
                 end if
                 if (somev) then
                    if (.not. select(ki)) cycle loop_260
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real left eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = ki + 1,n
                       work(k + iv*n) = -t(ki,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+1:n,ki+1:n) - wr ]**t * x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_ddot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          ! solve [ t(j,j) - wr ]**t * x = work
                          call la_dlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_dscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          vmax = max(abs(work(j + iv*n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_ddot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          work(j + 1 + iv*n) = work(j + 1 + iv*n) - la_ddot(j - ki - 1,t(ki + 1,j + 1) &
                                    ,1,work(ki + 1 + iv*n),1)
                          ! solve
                          ! [ t(j,j)-wr   t(j,j+1)      ]**t * x = scale*( work1 )
                          ! [ t(j+1,j)    t(j+1,j+1)-wr ]                ( work2 )
                          call la_dlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_dscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          work(j + 1 + iv*n) = x(2,1)
                          vmax = max(abs(work(j + iv*n)),abs(work(j + 1 + iv*n)),vmax)

                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_dcopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                       ii = la_idamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_dscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n) call la_dgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + iv*n),1,work(ki + iv*n),vl(1,ki),1)
                       ii = la_idamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_dscal(n,remax,vl(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex left eigenvector.
                    ! initial solve:
                    ! [ ( t(ki,ki)    t(ki,ki+1)  )**t - (wr - i* wi) ]*x = 0.
                    ! [ ( t(ki+1,ki) t(ki+1,ki+1) )                   ]
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + (iv)*n) = wi/t(ki,ki + 1)
                       work(ki + 1 + (iv + 1)*n) = one
                    else
                       work(ki + (iv)*n) = one
                       work(ki + 1 + (iv + 1)*n) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + (iv)*n) = zero
                    work(ki + (iv + 1)*n) = zero
                    ! form right-hand side.
                    do k = ki + 2,n
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(ki,k)
                       work(k + (iv + 1)*n) = -work(ki + 1 + (iv + 1)*n)*t(ki + 1,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+2:n,ki+2:n)**t - (wr-i*wi) ]*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_dscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_ddot(j - ki - 2,t(ki + 2,j) &
                                    ,1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_ddot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve [ t(j,j)-(wr-i*wi) ]*(x11+i*x12)= wk+i*wk2
                          call la_dlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_dscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          vmax = max(abs(work(j + (iv)*n)),abs(work(j + (iv + 1)*n)),vmax)

                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_dscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_dscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_ddot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_ddot(j - ki - 2,t(ki + 2, &
                                     j),1,work(ki + 2 + (iv + 1)*n),1)
                          work(j + 1 + (iv)*n) = work(j + 1 + (iv)*n) - la_ddot(j - ki - 2,t(ki + 2, &
                                     j + 1),1,work(ki + 2 + (iv)*n),1)
                          work(j + 1 + (iv + 1)*n) = work(j + 1 + (iv + 1)*n) - la_ddot(j - ki - 2,t(ki + &
                                    2,j + 1),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve 2-by-2 complex linear equation
                          ! [ (t(j,j)   t(j,j+1)  )**t - (wr-i*wi)*i ]*x = scale*b
                          ! [ (t(j+1,j) t(j+1,j+1))                  ]
                          call la_dlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_dscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_dscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          work(j + 1 + (iv)*n) = x(2,1)
                          work(j + 1 + (iv + 1)*n) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_dcopy(n - ki + 1,work(ki + (iv)*n),1,vl(ki,is),1)

                       call la_dcopy(n - ki + 1,work(ki + (iv + 1)*n),1,vl(ki,is + 1),1)

                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_dscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_dscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n - 1) then
                          call la_dgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv)*n),1,work(ki + (iv)*n),vl(1,ki),1)
                          call la_dgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv + 1)*n),1,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       else
                          call la_dscal(n,work(ki + (iv)*n),vl(1,ki),1)
                          call la_dscal(n,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_dscal(n,remax,vl(1,ki),1)
                       call la_dscal(n,remax,vl(1,ki + 1),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + (iv)*n) = zero
                          work(k + (iv + 1)*n) = zero
                       end do
                       iscomplex(iv) = ip
                       iscomplex(iv + 1) = -ip
                       iv = iv + 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki and ki+1)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki + 1
                    end if
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv >= nb - 1) .or. (ki2 == n)) then
                       call la_dgemm('N','N',n,iv,n - ki2 + iv,one,vl(1,ki2 - iv + 1),ldvl, &
                                 work(ki2 - iv + 1 + (1)*n),n,zero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_idamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_dscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_dlacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki2 - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if ! blocked back-transform
                 is = is + 1
                 if (ip /= 0) is = is + 1
              end do loop_260
           end if
           return
     end subroutine la_dtrevc3
#ifdef LA_WITH_XDP
     !> XTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by XHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**T)*T = w*(y**T)
     !> where y**T denotes the transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_xtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(xdp),intent(in) :: t(ldt,*)
           real(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,iv,maxwrk,nb, &
                     ki2
           real(xdp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(xdp) :: x(2,2)
           integer(ilp) :: iscomplex(nbmax)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           nb = la_ilaenv(1,'XTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           lquery = (lwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -14
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_xlaset('F',n,1 + 2*nb,zero,zero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_xlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_xlabad(unfl,ovfl)
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first  of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
             ! iscomplex array stores ip for each column in current block.
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! for complex right vector, uses iv-1 for real part and iv for complex part.
              ! non-blocked version always uses iv=2;
              ! blocked     version starts with iv=nb, goes down to 1 or 2.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 2
              if (nb > 2) then
                 iv = nb
              end if
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == -1) then
                    ! previous iteration (ki+1) was second of conjugate pair,
                    ! so this ki is first of conjugate pair; skip to end of loop
                    ip = 1
                    cycle loop_140
                 else if (ki == 1) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki,ki - 1) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is second of conjugate pair
                    ip = -1
                 end if
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) cycle loop_140
                    else
                       if (.not. select(ki - 1)) cycle loop_140
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real right eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = 1,ki - 1
                       work(k + iv*n) = -t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-1,1:ki-1) - wr ]*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_xlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_xscal(ki,scale,work(1 + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          ! update right-hand side
                          call la_xaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + iv*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_xlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_xscal(ki,scale,work(1 + iv*n),1)

                          work(j - 1 + iv*n) = x(1,1)
                          work(j + iv*n) = x(2,1)
                          ! update right-hand side
                          call la_xaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + iv*n),1)

                          call la_xaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + iv*n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_xcopy(ki,work(1 + iv*n),1,vr(1,is),1)
                       ii = la_ixamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_xscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 1) call la_xgemv('N',n,ki - 1,one,vr,ldvr,work(1 + iv*n), &
                                 1,work(ki + iv*n),vr(1,ki),1)
                       ii = la_ixamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_xscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex right eigenvector.
                    ! initial solve
                    ! [ ( t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i*wi) ]*x = 0.
                    ! [ ( t(ki,  ki-1) t(ki,  ki) )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + (iv - 1)*n) = one
                       work(ki + (iv)*n) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + (iv - 1)*n) = -wi/t(ki,ki - 1)
                       work(ki + (iv)*n) = one
                    end if
                    work(ki + (iv - 1)*n) = zero
                    work(ki - 1 + (iv)*n) = zero
                    ! form right-hand side.
                    do k = 1,ki - 2
                       work(k + (iv - 1)*n) = -work(ki - 1 + (iv - 1)*n)*t(k,ki - 1)
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-2,1:ki-2) - (wr+i*wi) ]*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_xlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_xscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j + (iv - 1)*n) = x(1,1)
                          work(j + (iv)*n) = x(1,2)
                          ! update the right-hand side
                          call la_xaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + (iv - 1)*n),1)

                          call la_xaxpy(j - 1,-x(1,2),t(1,j),1,work(1 + (iv)*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_xlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_xscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j - 1 + (iv - 1)*n) = x(1,1)
                          work(j + (iv - 1)*n) = x(2,1)
                          work(j - 1 + (iv)*n) = x(1,2)
                          work(j + (iv)*n) = x(2,2)
                          ! update the right-hand side
                          call la_xaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + (iv - 1)*n), &
                                     1)
                          call la_xaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + (iv - 1)*n), &
                                    1)
                          call la_xaxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + (iv)*n), &
                                    1)
                          call la_xaxpy(j - 2,-x(2,2),t(1,j),1,work(1 + (iv)*n),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_xcopy(ki,work(1 + (iv - 1)*n),1,vr(1,is - 1),1)
                       call la_xcopy(ki,work(1 + (iv)*n),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_xscal(ki,remax,vr(1,is - 1),1)
                       call la_xscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 2) then
                          call la_xgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv - 1)*n), &
                                    1,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_xgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv)*n),1, &
                                    work(ki + (iv)*n),vr(1,ki),1)
                       else
                          call la_xscal(n,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_xscal(n,work(ki + (iv)*n),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_xscal(n,remax,vr(1,ki - 1),1)
                       call la_xscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + (iv - 1)*n) = zero
                          work(k + (iv)*n) = zero
                       end do
                       iscomplex(iv - 1) = -ip
                       iscomplex(iv) = ip
                       iv = iv - 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki-1 and ki)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki - 1
                    end if
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv <= 2) .or. (ki2 == 1)) then
                       call la_xgemm('N','N',n,nb - iv + 1,ki2 + nb - iv,one,vr,ldvr,work(1 + &
                                 (iv)*n),n,zero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_ixamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_xscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_xlacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki2), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if ! blocked back-transform
                 is = is - 1
                 if (ip /= 0) is = is - 1
              end do loop_140
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! for complex left vector, uses iv for real part and iv+1 for complex part.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb-1 or nb.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 1
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == 1) then
                    ! previous iteration (ki-1) was first of conjugate pair,
                    ! so this ki is second of conjugate pair; skip to end of loop
                    ip = -1
                    cycle loop_260
                 else if (ki == n) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki + 1,ki) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is first of conjugate pair
                    ip = 1
                 end if
                 if (somev) then
                    if (.not. select(ki)) cycle loop_260
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real left eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = ki + 1,n
                       work(k + iv*n) = -t(ki,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+1:n,ki+1:n) - wr ]**t * x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_xdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          ! solve [ t(j,j) - wr ]**t * x = work
                          call la_xlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_xscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          vmax = max(abs(work(j + iv*n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_xdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          work(j + 1 + iv*n) = work(j + 1 + iv*n) - la_xdot(j - ki - 1,t(ki + 1,j + 1) &
                                    ,1,work(ki + 1 + iv*n),1)
                          ! solve
                          ! [ t(j,j)-wr   t(j,j+1)      ]**t * x = scale*( work1 )
                          ! [ t(j+1,j)    t(j+1,j+1)-wr ]                ( work2 )
                          call la_xlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_xscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          work(j + 1 + iv*n) = x(2,1)
                          vmax = max(abs(work(j + iv*n)),abs(work(j + 1 + iv*n)),vmax)

                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_xcopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                       ii = la_ixamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_xscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n) call la_xgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + iv*n),1,work(ki + iv*n),vl(1,ki),1)
                       ii = la_ixamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_xscal(n,remax,vl(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex left eigenvector.
                    ! initial solve:
                    ! [ ( t(ki,ki)    t(ki,ki+1)  )**t - (wr - i* wi) ]*x = 0.
                    ! [ ( t(ki+1,ki) t(ki+1,ki+1) )                   ]
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + (iv)*n) = wi/t(ki,ki + 1)
                       work(ki + 1 + (iv + 1)*n) = one
                    else
                       work(ki + (iv)*n) = one
                       work(ki + 1 + (iv + 1)*n) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + (iv)*n) = zero
                    work(ki + (iv + 1)*n) = zero
                    ! form right-hand side.
                    do k = ki + 2,n
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(ki,k)
                       work(k + (iv + 1)*n) = -work(ki + 1 + (iv + 1)*n)*t(ki + 1,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+2:n,ki+2:n)**t - (wr-i*wi) ]*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_xscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_xdot(j - ki - 2,t(ki + 2,j) &
                                    ,1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_xdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve [ t(j,j)-(wr-i*wi) ]*(x11+i*x12)= wk+i*wk2
                          call la_xlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_xscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          vmax = max(abs(work(j + (iv)*n)),abs(work(j + (iv + 1)*n)),vmax)

                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_xscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_xscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_xdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_xdot(j - ki - 2,t(ki + 2, &
                                     j),1,work(ki + 2 + (iv + 1)*n),1)
                          work(j + 1 + (iv)*n) = work(j + 1 + (iv)*n) - la_xdot(j - ki - 2,t(ki + 2, &
                                     j + 1),1,work(ki + 2 + (iv)*n),1)
                          work(j + 1 + (iv + 1)*n) = work(j + 1 + (iv + 1)*n) - la_xdot(j - ki - 2,t(ki + &
                                    2,j + 1),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve 2-by-2 complex linear equation
                          ! [ (t(j,j)   t(j,j+1)  )**t - (wr-i*wi)*i ]*x = scale*b
                          ! [ (t(j+1,j) t(j+1,j+1))                  ]
                          call la_xlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_xscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_xscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          work(j + 1 + (iv)*n) = x(2,1)
                          work(j + 1 + (iv + 1)*n) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_xcopy(n - ki + 1,work(ki + (iv)*n),1,vl(ki,is),1)

                       call la_xcopy(n - ki + 1,work(ki + (iv + 1)*n),1,vl(ki,is + 1),1)

                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_xscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_xscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n - 1) then
                          call la_xgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv)*n),1,work(ki + (iv)*n),vl(1,ki),1)
                          call la_xgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv + 1)*n),1,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       else
                          call la_xscal(n,work(ki + (iv)*n),vl(1,ki),1)
                          call la_xscal(n,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_xscal(n,remax,vl(1,ki),1)
                       call la_xscal(n,remax,vl(1,ki + 1),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + (iv)*n) = zero
                          work(k + (iv + 1)*n) = zero
                       end do
                       iscomplex(iv) = ip
                       iscomplex(iv + 1) = -ip
                       iv = iv + 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki and ki+1)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki + 1
                    end if
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv >= nb - 1) .or. (ki2 == n)) then
                       call la_xgemm('N','N',n,iv,n - ki2 + iv,one,vl(1,ki2 - iv + 1),ldvl, &
                                 work(ki2 - iv + 1 + (1)*n),n,zero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_ixamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_xscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_xlacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki2 - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if ! blocked back-transform
                 is = is + 1
                 if (ip /= 0) is = is + 1
              end do loop_260
           end if
           return
     end subroutine la_xtrevc3
#endif
#ifdef LA_WITH_QP
     !> QTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a real upper quasi-triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a real general matrix:  A = Q*T*Q**T, as computed by QHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**T)*T = w*(y**T)
     !> where y**T denotes the transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal blocks of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the orthogonal factor that reduces a matrix
     !> A to Schur form T, then Q*X and Q*Y are the matrices of right and
     !> left eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_qtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           real(qp),intent(in) :: t(ldt,*)
           real(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,pair,rightv,somev
           integer(ilp) :: i,ierr,ii,ip,is,j,j1,j2,jnxt,k,ki,iv,maxwrk,nb, &
                     ki2
           real(qp) :: beta,bignum,emax,ovfl,rec,remax,scale,smin,smlnum,ulp,unfl, &
                     vcrit,vmax,wi,wr,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Local Arrays
           real(qp) :: x(2,2)
           integer(ilp) :: iscomplex(nbmax)
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           info = 0
           nb = la_ilaenv(1,'QTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           lquery = (lwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -14
           else
              ! set m to the number of columns required to store the selected
              ! eigenvectors, standardize the array select if necessary, and
              ! test mm.
              if (somev) then
                 m = 0
                 pair = .false.
                 do j = 1,n
                    if (pair) then
                       pair = .false.
                       select(j) = .false.
                    else
                       if (j < n) then
                          if (t(j + 1,j) == zero) then
                             if (select(j)) m = m + 1
                          else
                             pair = .true.
                             if (select(j) .or. select(j + 1)) then
                                select(j) = .true.
                                m = m + 2
                             end if
                          end if
                       else
                          if (select(n)) m = m + 1
                       end if
                    end if
                 end do
              else
                 m = n
              end if
              if (mm < m) then
                 info = -11
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_qlaset('F',n,1 + 2*nb,zero,zero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_qlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_qlabad(unfl,ovfl)
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           work(1) = zero
           do j = 2,n
              work(j) = zero
              do i = 1,j - 1
                 work(j) = work(j) + abs(t(i,j))
              end do
           end do
           ! index ip is used to specify the real or complex eigenvalue:
             ! ip = 0, real eigenvalue,
                  ! 1, first  of conjugate complex pair: (wr,wi)
                 ! -1, second of conjugate complex pair: (wr,wi)
             ! iscomplex array stores ip for each column in current block.
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! for complex right vector, uses iv-1 for real part and iv for complex part.
              ! non-blocked version always uses iv=2;
              ! blocked     version starts with iv=nb, goes down to 1 or 2.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 2
              if (nb > 2) then
                 iv = nb
              end if
              ip = 0
              is = m
              loop_140: do ki = n,1,-1
                 if (ip == -1) then
                    ! previous iteration (ki+1) was second of conjugate pair,
                    ! so this ki is first of conjugate pair; skip to end of loop
                    ip = 1
                    cycle loop_140
                 else if (ki == 1) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki,ki - 1) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is second of conjugate pair
                    ip = -1
                 end if
                 if (somev) then
                    if (ip == 0) then
                       if (.not. select(ki)) cycle loop_140
                    else
                       if (.not. select(ki - 1)) cycle loop_140
                    end if
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki - 1)))*sqrt(abs(t(ki - 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real right eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = 1,ki - 1
                       work(k + iv*n) = -t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-1,1:ki-1) - wr ]*x = scale*work.
                    jnxt = ki - 1
                    loop_60: do j = ki - 1,1,-1
                       if (j > jnxt) cycle loop_60
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_qlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_qscal(ki,scale,work(1 + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          ! update right-hand side
                          call la_qaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + iv*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_qlaln2(.false.,2,1,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(2,1) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(2,1) = x(2,1)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) call la_qscal(ki,scale,work(1 + iv*n),1)

                          work(j - 1 + iv*n) = x(1,1)
                          work(j + iv*n) = x(2,1)
                          ! update right-hand side
                          call la_qaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + iv*n),1)

                          call la_qaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + iv*n),1)

                       end if
                    end do loop_60
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_qcopy(ki,work(1 + iv*n),1,vr(1,is),1)
                       ii = la_iqamax(ki,vr(1,is),1)
                       remax = one/abs(vr(ii,is))
                       call la_qscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 1) call la_qgemv('N',n,ki - 1,one,vr,ldvr,work(1 + iv*n), &
                                 1,work(ki + iv*n),vr(1,ki),1)
                       ii = la_iqamax(n,vr(1,ki),1)
                       remax = one/abs(vr(ii,ki))
                       call la_qscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex right eigenvector.
                    ! initial solve
                    ! [ ( t(ki-1,ki-1) t(ki-1,ki) ) - (wr + i*wi) ]*x = 0.
                    ! [ ( t(ki,  ki-1) t(ki,  ki) )               ]
                    if (abs(t(ki - 1,ki)) >= abs(t(ki,ki - 1))) then
                       work(ki - 1 + (iv - 1)*n) = one
                       work(ki + (iv)*n) = wi/t(ki - 1,ki)
                    else
                       work(ki - 1 + (iv - 1)*n) = -wi/t(ki,ki - 1)
                       work(ki + (iv)*n) = one
                    end if
                    work(ki + (iv - 1)*n) = zero
                    work(ki - 1 + (iv)*n) = zero
                    ! form right-hand side.
                    do k = 1,ki - 2
                       work(k + (iv - 1)*n) = -work(ki - 1 + (iv - 1)*n)*t(k,ki - 1)
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(k,ki)
                    end do
                    ! solve upper quasi-triangular system:
                    ! [ t(1:ki-2,1:ki-2) - (wr+i*wi) ]*x = scale*(work+i*work2)
                    jnxt = ki - 2
                    loop_90: do j = ki - 2,1,-1
                       if (j > jnxt) cycle loop_90
                       j1 = j
                       j2 = j
                       jnxt = j - 1
                       if (j > 1) then
                          if (t(j,j - 1) /= zero) then
                             j1 = j - 1
                             jnxt = j - 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          call la_qlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x(1,1) and x(1,2) to avoid overflow when
                          ! updating the right-hand side.
                          if (xnorm > one) then
                             if (work(j) > bignum/xnorm) then
                                x(1,1) = x(1,1)/xnorm
                                x(1,2) = x(1,2)/xnorm
                                scale = scale/xnorm
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_qscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j + (iv - 1)*n) = x(1,1)
                          work(j + (iv)*n) = x(1,2)
                          ! update the right-hand side
                          call la_qaxpy(j - 1,-x(1,1),t(1,j),1,work(1 + (iv - 1)*n),1)

                          call la_qaxpy(j - 1,-x(1,2),t(1,j),1,work(1 + (iv)*n),1)

                       else
                          ! 2-by-2 diagonal block
                          call la_qlaln2(.false.,2,2,smin,one,t(j - 1,j - 1),ldt,one, &
                                    one,work(j - 1 + (iv - 1)*n),n,wr,wi,x,2,scale,xnorm,ierr)
                          ! scale x to avoid overflow when updating
                          ! the right-hand side.
                          if (xnorm > one) then
                             beta = max(work(j - 1),work(j))
                             if (beta > bignum/xnorm) then
                                rec = one/xnorm
                                x(1,1) = x(1,1)*rec
                                x(1,2) = x(1,2)*rec
                                x(2,1) = x(2,1)*rec
                                x(2,2) = x(2,2)*rec
                                scale = scale*rec
                             end if
                          end if
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(ki,scale,work(1 + (iv - 1)*n),1)
                             call la_qscal(ki,scale,work(1 + (iv)*n),1)
                          end if
                          work(j - 1 + (iv - 1)*n) = x(1,1)
                          work(j + (iv - 1)*n) = x(2,1)
                          work(j - 1 + (iv)*n) = x(1,2)
                          work(j + (iv)*n) = x(2,2)
                          ! update the right-hand side
                          call la_qaxpy(j - 2,-x(1,1),t(1,j - 1),1,work(1 + (iv - 1)*n), &
                                     1)
                          call la_qaxpy(j - 2,-x(2,1),t(1,j),1,work(1 + (iv - 1)*n), &
                                    1)
                          call la_qaxpy(j - 2,-x(1,2),t(1,j - 1),1,work(1 + (iv)*n), &
                                    1)
                          call la_qaxpy(j - 2,-x(2,2),t(1,j),1,work(1 + (iv)*n),1)

                       end if
                    end do loop_90
                    ! copy the vector x or q*x to vr and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vr and normalize.
                       call la_qcopy(ki,work(1 + (iv - 1)*n),1,vr(1,is - 1),1)
                       call la_qcopy(ki,work(1 + (iv)*n),1,vr(1,is),1)
                       emax = zero
                       do k = 1,ki
                          emax = max(emax,abs(vr(k,is - 1)) + abs(vr(k,is)))
                       end do
                       remax = one/emax
                       call la_qscal(ki,remax,vr(1,is - 1),1)
                       call la_qscal(ki,remax,vr(1,is),1)
                       do k = ki + 1,n
                          vr(k,is - 1) = zero
                          vr(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki > 2) then
                          call la_qgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv - 1)*n), &
                                    1,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_qgemv('N',n,ki - 2,one,vr,ldvr,work(1 + (iv)*n),1, &
                                    work(ki + (iv)*n),vr(1,ki),1)
                       else
                          call la_qscal(n,work(ki - 1 + (iv - 1)*n),vr(1,ki - 1),1)
                          call la_qscal(n,work(ki + (iv)*n),vr(1,ki),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vr(k,ki - 1)) + abs(vr(k,ki)))
                       end do
                       remax = one/emax
                       call la_qscal(n,remax,vr(1,ki - 1),1)
                       call la_qscal(n,remax,vr(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out below vector
                       do k = ki + 1,n
                          work(k + (iv - 1)*n) = zero
                          work(k + (iv)*n) = zero
                       end do
                       iscomplex(iv - 1) = -ip
                       iscomplex(iv) = ip
                       iv = iv - 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki-1 and ki)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki - 1
                    end if
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv <= 2) .or. (ki2 == 1)) then
                       call la_qgemm('N','N',n,nb - iv + 1,ki2 + nb - iv,one,vr,ldvr,work(1 + &
                                 (iv)*n),n,zero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_iqamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_qscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_qlacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki2), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if ! blocked back-transform
                 is = is - 1
                 if (ip /= 0) is = is - 1
              end do loop_140
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! for complex left vector, uses iv for real part and iv+1 for complex part.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb-1 or nb.
              ! (note the "0-th" column is used for 1-norms computed above.)
              iv = 1
              ip = 0
              is = 1
              loop_260: do ki = 1,n
                 if (ip == 1) then
                    ! previous iteration (ki-1) was first of conjugate pair,
                    ! so this ki is second of conjugate pair; skip to end of loop
                    ip = -1
                    cycle loop_260
                 else if (ki == n) then
                    ! last column, so this ki must be real eigenvalue
                    ip = 0
                 else if (t(ki + 1,ki) == zero) then
                    ! zero on sub-diagonal, so this ki is real eigenvalue
                    ip = 0
                 else
                    ! non-zero on sub-diagonal, so this ki is first of conjugate pair
                    ip = 1
                 end if
                 if (somev) then
                    if (.not. select(ki)) cycle loop_260
                 end if
                 ! compute the ki-th eigenvalue (wr,wi).
                 wr = t(ki,ki)
                 wi = zero
                 if (ip /= 0) wi = sqrt(abs(t(ki,ki + 1)))*sqrt(abs(t(ki + 1,ki)))
                 smin = max(ulp*(abs(wr) + abs(wi)),smlnum)
                 if (ip == 0) then
                    ! --------------------------------------------------------
                    ! real left eigenvector
                    work(ki + iv*n) = one
                    ! form right-hand side.
                    do k = ki + 1,n
                       work(k + iv*n) = -t(ki,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+1:n,ki+1:n) - wr ]**t * x = scale*work
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 1
                    loop_170: do j = ki + 1,n
                       if (j < jnxt) cycle loop_170
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_qdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          ! solve [ t(j,j) - wr ]**t * x = work
                          call la_qlaln2(.false.,1,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_qscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          vmax = max(abs(work(j + iv*n)),vmax)
                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + iv*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + iv*n) = work(j + iv*n) - la_qdot(j - ki - 1,t(ki + 1,j),1, &
                                    work(ki + 1 + iv*n),1)
                          work(j + 1 + iv*n) = work(j + 1 + iv*n) - la_qdot(j - ki - 1,t(ki + 1,j + 1) &
                                    ,1,work(ki + 1 + iv*n),1)
                          ! solve
                          ! [ t(j,j)-wr   t(j,j+1)      ]**t * x = scale*( work1 )
                          ! [ t(j+1,j)    t(j+1,j+1)-wr ]                ( work2 )
                          call la_qlaln2(.true.,2,1,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,zero,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) call la_qscal(n - ki + 1,scale,work(ki + iv*n),1)

                          work(j + iv*n) = x(1,1)
                          work(j + 1 + iv*n) = x(2,1)
                          vmax = max(abs(work(j + iv*n)),abs(work(j + 1 + iv*n)),vmax)

                          vcrit = bignum/vmax
                       end if
                    end do loop_170
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_qcopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                       ii = la_iqamax(n - ki + 1,vl(ki,is),1) + ki - 1
                       remax = one/abs(vl(ii,is))
                       call la_qscal(n - ki + 1,remax,vl(ki,is),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n) call la_qgemv('N',n,n - ki,one,vl(1,ki + 1),ldvl,work( &
                                 ki + 1 + iv*n),1,work(ki + iv*n),vl(1,ki),1)
                       ii = la_iqamax(n,vl(1,ki),1)
                       remax = one/abs(vl(ii,ki))
                       call la_qscal(n,remax,vl(1,ki),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + iv*n) = zero
                       end do
                       iscomplex(iv) = ip
                       ! back-transform and normalization is done below
                    end if
                 else
                    ! --------------------------------------------------------
                    ! complex left eigenvector.
                    ! initial solve:
                    ! [ ( t(ki,ki)    t(ki,ki+1)  )**t - (wr - i* wi) ]*x = 0.
                    ! [ ( t(ki+1,ki) t(ki+1,ki+1) )                   ]
                    if (abs(t(ki,ki + 1)) >= abs(t(ki + 1,ki))) then
                       work(ki + (iv)*n) = wi/t(ki,ki + 1)
                       work(ki + 1 + (iv + 1)*n) = one
                    else
                       work(ki + (iv)*n) = one
                       work(ki + 1 + (iv + 1)*n) = -wi/t(ki + 1,ki)
                    end if
                    work(ki + 1 + (iv)*n) = zero
                    work(ki + (iv + 1)*n) = zero
                    ! form right-hand side.
                    do k = ki + 2,n
                       work(k + (iv)*n) = -work(ki + (iv)*n)*t(ki,k)
                       work(k + (iv + 1)*n) = -work(ki + 1 + (iv + 1)*n)*t(ki + 1,k)
                    end do
                    ! solve transposed quasi-triangular system:
                    ! [ t(ki+2:n,ki+2:n)**t - (wr-i*wi) ]*x = work1+i*work2
                    vmax = one
                    vcrit = bignum
                    jnxt = ki + 2
                    loop_200: do j = ki + 2,n
                       if (j < jnxt) cycle loop_200
                       j1 = j
                       j2 = j
                       jnxt = j + 1
                       if (j < n) then
                          if (t(j + 1,j) /= zero) then
                             j2 = j + 1
                             jnxt = j + 2
                          end if
                       end if
                       if (j1 == j2) then
                          ! 1-by-1 diagonal block
                          ! scale if necessary to avoid overflow when
                          ! forming the right-hand side elements.
                          if (work(j) > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_qscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_qdot(j - ki - 2,t(ki + 2,j) &
                                    ,1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_qdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve [ t(j,j)-(wr-i*wi) ]*(x11+i*x12)= wk+i*wk2
                          call la_qlaln2(.false.,1,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_qscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          vmax = max(abs(work(j + (iv)*n)),abs(work(j + (iv + 1)*n)),vmax)

                          vcrit = bignum/vmax
                       else
                          ! 2-by-2 diagonal block
                          ! scale if necessary to avoid overflow when forming
                          ! the right-hand side elements.
                          beta = max(work(j),work(j + 1))
                          if (beta > vcrit) then
                             rec = one/vmax
                             call la_qscal(n - ki + 1,rec,work(ki + (iv)*n),1)
                             call la_qscal(n - ki + 1,rec,work(ki + (iv + 1)*n),1)
                             vmax = one
                             vcrit = bignum
                          end if
                          work(j + (iv)*n) = work(j + (iv)*n) - la_qdot(j - ki - 2,t(ki + 2, &
                                    j),1,work(ki + 2 + (iv)*n),1)
                          work(j + (iv + 1)*n) = work(j + (iv + 1)*n) - la_qdot(j - ki - 2,t(ki + 2, &
                                     j),1,work(ki + 2 + (iv + 1)*n),1)
                          work(j + 1 + (iv)*n) = work(j + 1 + (iv)*n) - la_qdot(j - ki - 2,t(ki + 2, &
                                     j + 1),1,work(ki + 2 + (iv)*n),1)
                          work(j + 1 + (iv + 1)*n) = work(j + 1 + (iv + 1)*n) - la_qdot(j - ki - 2,t(ki + &
                                    2,j + 1),1,work(ki + 2 + (iv + 1)*n),1)
                          ! solve 2-by-2 complex linear equation
                          ! [ (t(j,j)   t(j,j+1)  )**t - (wr-i*wi)*i ]*x = scale*b
                          ! [ (t(j+1,j) t(j+1,j+1))                  ]
                          call la_qlaln2(.true.,2,2,smin,one,t(j,j),ldt,one,one, &
                                    work(j + iv*n),n,wr,-wi,x,2,scale,xnorm,ierr)
                          ! scale if necessary
                          if (scale /= one) then
                             call la_qscal(n - ki + 1,scale,work(ki + (iv)*n),1)
                             call la_qscal(n - ki + 1,scale,work(ki + (iv + 1)*n),1)
                          end if
                          work(j + (iv)*n) = x(1,1)
                          work(j + (iv + 1)*n) = x(1,2)
                          work(j + 1 + (iv)*n) = x(2,1)
                          work(j + 1 + (iv + 1)*n) = x(2,2)
                          vmax = max(abs(x(1,1)),abs(x(1,2)),abs(x(2,1)),abs(x( &
                                     2,2)),vmax)
                          vcrit = bignum/vmax
                       end if
                    end do loop_200
                    ! copy the vector x or q*x to vl and normalize.
                    if (.not. over) then
                       ! ------------------------------
                       ! no back-transform: copy x to vl and normalize.
                       call la_qcopy(n - ki + 1,work(ki + (iv)*n),1,vl(ki,is),1)

                       call la_qcopy(n - ki + 1,work(ki + (iv + 1)*n),1,vl(ki,is + 1),1)

                       emax = zero
                       do k = ki,n
                          emax = max(emax,abs(vl(k,is)) + abs(vl(k,is + 1)))
                       end do
                       remax = one/emax
                       call la_qscal(n - ki + 1,remax,vl(ki,is),1)
                       call la_qscal(n - ki + 1,remax,vl(ki,is + 1),1)
                       do k = 1,ki - 1
                          vl(k,is) = zero
                          vl(k,is + 1) = zero
                       end do
                    else if (nb == 1) then
                       ! ------------------------------
                       ! version 1: back-transform each vector with gemv, q*x.
                       if (ki < n - 1) then
                          call la_qgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv)*n),1,work(ki + (iv)*n),vl(1,ki),1)
                          call la_qgemv('N',n,n - ki - 1,one,vl(1,ki + 2),ldvl,work(ki + 2 + &
                                    (iv + 1)*n),1,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       else
                          call la_qscal(n,work(ki + (iv)*n),vl(1,ki),1)
                          call la_qscal(n,work(ki + 1 + (iv + 1)*n),vl(1,ki + 1),1)
                       end if
                       emax = zero
                       do k = 1,n
                          emax = max(emax,abs(vl(k,ki)) + abs(vl(k,ki + 1)))
                       end do
                       remax = one/emax
                       call la_qscal(n,remax,vl(1,ki),1)
                       call la_qscal(n,remax,vl(1,ki + 1),1)
                    else
                       ! ------------------------------
                       ! version 2: back-transform block of vectors with gemm
                       ! zero out above vector
                       ! could go from ki-nv+1 to ki-1
                       do k = 1,ki - 1
                          work(k + (iv)*n) = zero
                          work(k + (iv + 1)*n) = zero
                       end do
                       iscomplex(iv) = ip
                       iscomplex(iv + 1) = -ip
                       iv = iv + 1
                       ! back-transform and normalization is done below
                    end if
                 end if
                 if (nb > 1) then
                    ! --------------------------------------------------------
                    ! blocked version of back-transform
                    ! for complex case, ki2 includes both vectors (ki and ki+1)
                    if (ip == 0) then
                       ki2 = ki
                    else
                       ki2 = ki + 1
                    end if
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb-1 or nb,
                    ! or if this was last vector, do the gemm
                    if ((iv >= nb - 1) .or. (ki2 == n)) then
                       call la_qgemm('N','N',n,iv,n - ki2 + iv,one,vl(1,ki2 - iv + 1),ldvl, &
                                 work(ki2 - iv + 1 + (1)*n),n,zero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          if (iscomplex(k) == 0) then
                             ! real eigenvector
                             ii = la_iqamax(n,work(1 + (nb + k)*n),1)
                             remax = one/abs(work(ii + (nb + k)*n))
                          else if (iscomplex(k) == 1) then
                             ! first eigenvector of conjugate pair
                             emax = zero
                             do ii = 1,n
                                emax = max(emax,abs(work(ii + (nb + k)*n)) + abs(work(ii + ( &
                                          nb + k + 1)*n)))
                             end do
                             remax = one/emax
                          ! else if iscomplex(k)==-1
                             ! second eigenvector of conjugate pair
                             ! reuse same remax as previous k
                          end if
                          call la_qscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_qlacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki2 - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if ! blocked back-transform
                 is = is + 1
                 if (ip /= 0) is = is + 1
              end do loop_260
           end if
           return
     end subroutine la_qtrevc3
#endif

     !> STRSYL: solves the real Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**T, and  A and B are both upper quasi-
     !> triangular. A is M-by-M and B is N-by-N; the right hand side C and
     !> the solution X are M-by-N; and scale is an output scale factor, set
     !> <= 1 to avoid overflow in X.
     !> A and B must be in Schur canonical form (as returned by SHSEQR), that
     !> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
     !> each 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_strsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_sp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: ierr,j,k,k1,k2,knext,l,l1,l2,lnext
           real(sp) :: a11,bignum,da11,db,eps,scaloc,sgn,smin,smlnum,suml,sumr, &
                     xnorm
           ! Local Arrays
           real(sp) :: dum(1),vec(2,2),x(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'T') .and. .not. la_lsame(trana, &
                     'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'T') .and. .not. la_lsame( &
                     tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('STRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=sp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_slange('M',m,m,a,lda,dum),eps*la_slange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
               ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                        ! m                         l-1
              ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)].
                      ! i=k+1                       j=1
              ! start column loop (index = l)
              ! l1 (l2) : column index of the first (first) row of x(k,l).
              lnext = 1
              loop_70: do l = 1,n
                 if (l < lnext) cycle loop_70
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l).
                 knext = m
                 loop_60: do k = m,1,-1
                    if (k > knext) cycle loop_60
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_slaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_slaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_slasy2(.false.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
              end do loop_60
              end do loop_70
           else if (.not. notrna .and. notrnb) then
              ! solve    a**t *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1                          l-1
                ! r(k,l) = sum [a(i,k)**t*x(i,l)] +isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                          j=1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = 1
              loop_130: do l = 1,n
                 if (l < lnext) cycle loop_130
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_120: do k = 1,m
                    if (k < knext) cycle loop_120
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_slaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_slaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_sdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_sdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_slasy2(.true.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
              end do loop_120
              end do loop_130
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**t*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! top-right corner column by column by
                 ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                           ! k-1                            n
                  ! r(k,l) = sum [a(i,k)**t*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                           ! i=1                          j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_190: do l = n,1,-1
                 if (l > lnext) cycle loop_190
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_180: do k = 1,m
                    if (k < knext) cycle loop_180
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_slaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_slaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_sdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n) &
                                  ),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_slasy2(.true.,.true.,isgn,2,2,a(k1,k1),lda,b(l1,l1 &
                                 ),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
              end do loop_180
              end do loop_190
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-right corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                            ! m                          n
                  ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                          ! i=k+1                      j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_250: do l = n,1,-1
                 if (l > lnext) cycle loop_250
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = m
                 loop_240: do k = m,1,-1
                    if (k > knext) cycle loop_240
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_slaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_sdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_slaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_sdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_sdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_slasy2(.false.,.true.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_sscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
              end do loop_240
              end do loop_250
           end if
           return
     end subroutine la_strsyl
     !> DTRSYL: solves the real Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**T, and  A and B are both upper quasi-
     !> triangular. A is M-by-M and B is N-by-N; the right hand side C and
     !> the solution X are M-by-N; and scale is an output scale factor, set
     !> <= 1 to avoid overflow in X.
     !> A and B must be in Schur canonical form (as returned by DHSEQR), that
     !> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
     !> each 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_dtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_dp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: ierr,j,k,k1,k2,knext,l,l1,l2,lnext
           real(dp) :: a11,bignum,da11,db,eps,scaloc,sgn,smin,smlnum,suml,sumr, &
                     xnorm
           ! Local Arrays
           real(dp) :: dum(1),vec(2,2),x(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'T') .and. .not. la_lsame(trana, &
                     'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'T') .and. .not. la_lsame( &
                     tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=dp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_dlange('M',m,m,a,lda,dum),eps*la_dlange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
               ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                        ! m                         l-1
              ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)].
                      ! i=k+1                       j=1
              ! start column loop (index = l)
              ! l1 (l2) : column index of the first (first) row of x(k,l).
              lnext = 1
              loop_60: do l = 1,n
                 if (l < lnext) cycle loop_60
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l).
                 knext = m
                 loop_50: do k = m,1,-1
                    if (k > knext) cycle loop_50
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_dlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_dlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_dlasy2(.false.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_50
              end do loop_60
           else if (.not. notrna .and. notrnb) then
              ! solve    a**t *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1        t                    l-1
                ! r(k,l) = sum [a(i,k)**t*x(i,l)] +isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                          j=1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = 1
              loop_120: do l = 1,n
                 if (l < lnext) cycle loop_120
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_110: do k = 1,m
                    if (k < knext) cycle loop_110
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_dlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_dlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_ddot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_ddot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_dlasy2(.true.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_110
              end do loop_120
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**t*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! top-right corner column by column by
                 ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                           ! k-1                            n
                  ! r(k,l) = sum [a(i,k)**t*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                           ! i=1                          j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_180: do l = n,1,-1
                 if (l > lnext) cycle loop_180
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_170: do k = 1,m
                    if (k < knext) cycle loop_170
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_dlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_dlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_ddot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_dlasy2(.true.,.true.,isgn,2,2,a(k1,k1),lda,b(l1,l1 &
                                 ),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_170
              end do loop_180
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-right corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                            ! m                          n
                  ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                          ! i=k+1                      j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_240: do l = n,1,-1
                 if (l > lnext) cycle loop_240
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = m
                 loop_230: do k = m,1,-1
                    if (k > knext) cycle loop_230
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_dlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_ddot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_dlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_ddot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_ddot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_dlasy2(.false.,.true.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_dscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_230
              end do loop_240
           end if
           return
     end subroutine la_dtrsyl
#ifdef LA_WITH_XDP
     !> XTRSYL: solves the real Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**T, and  A and B are both upper quasi-
     !> triangular. A is M-by-M and B is N-by-N; the right hand side C and
     !> the solution X are M-by-N; and scale is an output scale factor, set
     !> <= 1 to avoid overflow in X.
     !> A and B must be in Schur canonical form (as returned by XHSEQR), that
     !> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
     !> each 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_xtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_xdp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(xdp),intent(out) :: scale
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: ierr,j,k,k1,k2,knext,l,l1,l2,lnext
           real(xdp) :: a11,bignum,da11,db,eps,scaloc,sgn,smin,smlnum,suml,sumr, &
                     xnorm
           ! Local Arrays
           real(xdp) :: dum(1),vec(2,2),x(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'T') .and. .not. la_lsame(trana, &
                     'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'T') .and. .not. la_lsame( &
                     tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('XTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=xdp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_xlange('M',m,m,a,lda,dum),eps*la_xlange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
               ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                        ! m                         l-1
              ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)].
                      ! i=k+1                       j=1
              ! start column loop (index = l)
              ! l1 (l2) : column index of the first (first) row of x(k,l).
              lnext = 1
              loop_60: do l = 1,n
                 if (l < lnext) cycle loop_60
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l).
                 knext = m
                 loop_50: do k = m,1,-1
                    if (k > knext) cycle loop_50
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_xlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_xlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_xlasy2(.false.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_50
              end do loop_60
           else if (.not. notrna .and. notrnb) then
              ! solve    a**t *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1        t                    l-1
                ! r(k,l) = sum [a(i,k)**t*x(i,l)] +isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                          j=1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = 1
              loop_120: do l = 1,n
                 if (l < lnext) cycle loop_120
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_110: do k = 1,m
                    if (k < knext) cycle loop_110
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_xlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_xlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_xdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_xdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_xlasy2(.true.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_110
              end do loop_120
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**t*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! top-right corner column by column by
                 ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                           ! k-1                            n
                  ! r(k,l) = sum [a(i,k)**t*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                           ! i=1                          j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_180: do l = n,1,-1
                 if (l > lnext) cycle loop_180
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_170: do k = 1,m
                    if (k < knext) cycle loop_170
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_xlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_xlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_xdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_xlasy2(.true.,.true.,isgn,2,2,a(k1,k1),lda,b(l1,l1 &
                                 ),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_170
              end do loop_180
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-right corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                            ! m                          n
                  ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                          ! i=k+1                      j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_240: do l = n,1,-1
                 if (l > lnext) cycle loop_240
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = m
                 loop_230: do k = m,1,-1
                    if (k > knext) cycle loop_230
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_xlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_xdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_xlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_xdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_xdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_xlasy2(.false.,.true.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_xscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_230
              end do loop_240
           end if
           return
     end subroutine la_xtrsyl
#endif
#ifdef LA_WITH_QP
     !> QTRSYL: solves the real Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**T, and  A and B are both upper quasi-
     !> triangular. A is M-by-M and B is N-by-N; the right hand side C and
     !> the solution X are M-by-N; and scale is an output scale factor, set
     !> <= 1 to avoid overflow in X.
     !> A and B must be in Schur canonical form (as returned by QHSEQR), that
     !> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
     !> each 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_qtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_qp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: ierr,j,k,k1,k2,knext,l,l1,l2,lnext
           real(qp) :: a11,bignum,da11,db,eps,scaloc,sgn,smin,smlnum,suml,sumr, &
                     xnorm
           ! Local Arrays
           real(qp) :: dum(1),vec(2,2),x(2,2)
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'T') .and. .not. la_lsame(trana, &
                     'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'T') .and. .not. la_lsame( &
                     tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=qp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_qlange('M',m,m,a,lda,dum),eps*la_qlange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
               ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                        ! m                         l-1
              ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)].
                      ! i=k+1                       j=1
              ! start column loop (index = l)
              ! l1 (l2) : column index of the first (first) row of x(k,l).
              lnext = 1
              loop_60: do l = 1,n
                 if (l < lnext) cycle loop_60
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l).
                 knext = m
                 loop_50: do k = m,1,-1
                    if (k > knext) cycle loop_50
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_qlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_qlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_qlasy2(.false.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_50
              end do loop_60
           else if (.not. notrna .and. notrnb) then
              ! solve    a**t *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1        t                    l-1
                ! r(k,l) = sum [a(i,k)**t*x(i,l)] +isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                          j=1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = 1
              loop_120: do l = 1,n
                 if (l < lnext) cycle loop_120
                 if (l == n) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l + 1,l) /= zero) then
                       l1 = l
                       l2 = l + 1
                       lnext = l + 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l + 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_110: do k = 1,m
                    if (k < knext) cycle loop_110
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_qlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_qlaln2(.true.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l1),1)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_qdot(l1 - 1,c(k1,1),ldc,b(1,l2),1)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l1),1)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_qdot(l1 - 1,c(k2,1),ldc,b(1,l2),1)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_qlasy2(.true.,.false.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_110
              end do loop_120
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**t*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! top-right corner column by column by
                 ! a(k,k)**t*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                           ! k-1                            n
                  ! r(k,l) = sum [a(i,k)**t*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                           ! i=1                          j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_180: do l = n,1,-1
                 if (l > lnext) cycle loop_180
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = 1
                 loop_170: do k = 1,m
                    if (k < knext) cycle loop_170
                    if (k == m) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k + 1,k) /= zero) then
                          k1 = k
                          k2 = k + 1
                          knext = k + 2
                       else
                          k1 = k
                          k2 = k
                          knext = k + 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_qlaln2(.true.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_qlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k1),1,c(1,l2),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l1),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_qdot(k1 - 1,a(1,k2),1,c(1,l2),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_qlasy2(.true.,.true.,isgn,2,2,a(k1,k1),lda,b(l1,l1 &
                                 ),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_170
              end do loop_180
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**t = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-right corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l)**t = c(k,l) - r(k,l)
              ! where
                            ! m                          n
                  ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b(l,j)**t].
                          ! i=k+1                      j=l+1
              ! start column loop (index = l)
              ! l1 (l2): column index of the first (last) row of x(k,l)
              lnext = n
              loop_240: do l = n,1,-1
                 if (l > lnext) cycle loop_240
                 if (l == 1) then
                    l1 = l
                    l2 = l
                 else
                    if (b(l,l - 1) /= zero) then
                       l1 = l - 1
                       l2 = l
                       lnext = l - 2
                    else
                       l1 = l
                       l2 = l
                       lnext = l - 1
                    end if
                 end if
                 ! start row loop (index = k)
                 ! k1 (k2): row index of the first (last) row of x(k,l)
                 knext = m
                 loop_230: do k = m,1,-1
                    if (k > knext) cycle loop_230
                    if (k == 1) then
                       k1 = k
                       k2 = k
                    else
                       if (a(k,k - 1) /= zero) then
                          k1 = k - 1
                          k2 = k
                          knext = k - 2
                       else
                          k1 = k
                          k2 = k
                          knext = k - 1
                       end if
                    end if
                    if (l1 == l2 .and. k1 == k2) then
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l1,c(k1,min(l1 + 1,n)),ldc,b(l1,min(l1 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       scaloc = one
                       a11 = a(k1,k1) + sgn*b(l1,l1)
                       da11 = abs(a11)
                       if (da11 <= smin) then
                          a11 = smin
                          da11 = smin
                          info = 1
                       end if
                       db = abs(vec(1,1))
                       if (da11 < one .and. db > one) then
                          if (db > bignum*da11) scaloc = one/db
                       end if
                       x(1,1) = (vec(1,1)*scaloc)/a11
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                    else if (l1 == l2 .and. k1 /= k2) then
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       call la_qlaln2(.false.,2,1,smin,one,a(k1,k1),lda,one,one, &
                                 vec,2,-sgn*b(l1,l1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k2,l1) = x(2,1)
                    else if (l1 /= l2 .and. k1 == k2) then
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = sgn*(c(k1,l1) - (suml + sgn*sumr))
                       suml = la_qdot(m - k1,a(k1,min(k1 + 1,m)),lda,c(min(k1 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = sgn*(c(k1,l2) - (suml + sgn*sumr))
                       call la_qlaln2(.false.,2,1,smin,one,b(l1,l1),ldb,one,one, &
                                 vec,2,-sgn*a(k1,k1),zero,x,2,scaloc,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(2,1)
                    else if (l1 /= l2 .and. k1 /= k2) then
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,1) = c(k1,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k1,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(n - l2,c(k1,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(1,2) = c(k1,l2) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l1),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l1,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,1) = c(k2,l1) - (suml + sgn*sumr)
                       suml = la_qdot(m - k2,a(k2,min(k2 + 1,m)),lda,c(min(k2 + 1,m), &
                                 l2),1)
                       sumr = la_qdot(n - l2,c(k2,min(l2 + 1,n)),ldc,b(l2,min(l2 + 1,n &
                                 )),ldb)
                       vec(2,2) = c(k2,l2) - (suml + sgn*sumr)
                       call la_qlasy2(.false.,.true.,isgn,2,2,a(k1,k1),lda,b(l1, &
                                 l1),ldb,vec,2,scaloc,x,2,xnorm,ierr)
                       if (ierr /= 0) info = 1
                       if (scaloc /= one) then
                          do j = 1,n
                             call la_qscal(m,scaloc,c(1,j),1)
                          end do
                          scale = scale*scaloc
                       end if
                       c(k1,l1) = x(1,1)
                       c(k1,l2) = x(1,2)
                       c(k2,l1) = x(2,1)
                       c(k2,l2) = x(2,2)
                    end if
                 end do loop_230
              end do loop_240
           end if
           return
     end subroutine la_qtrsyl
#endif

     !> SHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a real upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_shsein(side,eigsrc,initv,select,n,h,ldh,wr,wi,vl,ldvl,vr,ldvr, &
               mm,m,work,ifaill,ifailr,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(sp),intent(in) :: h(ldh,*),wi(*)
           real(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),wr(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,pair,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ksi,ksr,ldwork
           real(sp) :: bignum,eps3,hnorm,smlnum,ulp,unfl,wki,wkr
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors, and standardize the array select.
           m = 0
           pair = .false.
           do k = 1,n
              if (pair) then
                 pair = .false.
                 select(k) = .false.
              else
                 if (wi(k) == zero) then
                    if (select(k)) m = m + 1
                 else
                    pair = .true.
                    if (select(k) .or. select(k + 1)) then
                       select(k) = .true.
                       m = m + 2
                    end if
                 end if
              end if
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -13
           else if (mm < m) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('SHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_slamch('SAFE MINIMUM')
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ldwork = n + 1
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ksr = 1
           loop_120: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are zero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == zero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == zero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_slanhs('I',kr - kl + 1,h(kl,kl),ldh,work)
                    if (la_sisnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > zero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wkr = wr(k)
                 wki = wi(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. abs(wr(i) - wkr) + abs(wi(i) - wki) < eps3) &
                              then
                       wkr = wkr + eps3
                       go to 60
                    end if
                 end do
                 wr(k) = wkr
                 pair = wki /= zero
                 if (pair) then
                    ksi = ksr + 1
                 else
                    ksi = ksr
                 end if
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_slaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wkr,wki,vl( &
                    kl,ksr),vl(kl,ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum, &
                              iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifaill(ksr) = k
                       ifaill(ksi) = k
                    else
                       ifaill(ksr) = 0
                       ifaill(ksi) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = 1,kl - 1
                          vl(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_slaein(.true.,noinit,kr,h,ldh,wkr,wki,vr(1,ksr),vr(1, &
                              ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum,iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifailr(ksr) = k
                       ifailr(ksi) = k
                    else
                       ifailr(ksr) = 0
                       ifailr(ksi) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = kr + 1,n
                          vr(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (pair) then
                    ksr = ksr + 2
                 else
                    ksr = ksr + 1
                 end if
              end if
           end do loop_120
           return
     end subroutine la_shsein
     !> DHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a real upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_dhsein(side,eigsrc,initv,select,n,h,ldh,wr,wi,vl,ldvl,vr,ldvr, &
               mm,m,work,ifaill,ifailr,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(dp),intent(in) :: h(ldh,*),wi(*)
           real(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),wr(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,pair,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ksi,ksr,ldwork
           real(dp) :: bignum,eps3,hnorm,smlnum,ulp,unfl,wki,wkr
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors, and standardize the array select.
           m = 0
           pair = .false.
           do k = 1,n
              if (pair) then
                 pair = .false.
                 select(k) = .false.
              else
                 if (wi(k) == zero) then
                    if (select(k)) m = m + 1
                 else
                    pair = .true.
                    if (select(k) .or. select(k + 1)) then
                       select(k) = .true.
                       m = m + 2
                    end if
                 end if
              end if
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -13
           else if (mm < m) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('DHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_dlamch('SAFE MINIMUM')
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ldwork = n + 1
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ksr = 1
           loop_120: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are zero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == zero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == zero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_dlanhs('I',kr - kl + 1,h(kl,kl),ldh,work)
                    if (la_disnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > zero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wkr = wr(k)
                 wki = wi(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. abs(wr(i) - wkr) + abs(wi(i) - wki) < eps3) &
                              then
                       wkr = wkr + eps3
                       go to 60
                    end if
                 end do
                 wr(k) = wkr
                 pair = wki /= zero
                 if (pair) then
                    ksi = ksr + 1
                 else
                    ksi = ksr
                 end if
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_dlaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wkr,wki,vl( &
                    kl,ksr),vl(kl,ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum, &
                              iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifaill(ksr) = k
                       ifaill(ksi) = k
                    else
                       ifaill(ksr) = 0
                       ifaill(ksi) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = 1,kl - 1
                          vl(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_dlaein(.true.,noinit,kr,h,ldh,wkr,wki,vr(1,ksr),vr(1, &
                              ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum,iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifailr(ksr) = k
                       ifailr(ksi) = k
                    else
                       ifailr(ksr) = 0
                       ifailr(ksi) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = kr + 1,n
                          vr(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (pair) then
                    ksr = ksr + 2
                 else
                    ksr = ksr + 1
                 end if
              end if
           end do loop_120
           return
     end subroutine la_dhsein
#ifdef LA_WITH_XDP
     !> XHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a real upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_xhsein(side,eigsrc,initv,select,n,h,ldh,wr,wi,vl,ldvl,vr,ldvr, &
               mm,m,work,ifaill,ifailr,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(xdp),intent(in) :: h(ldh,*),wi(*)
           real(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),wr(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,pair,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ksi,ksr,ldwork
           real(xdp) :: bignum,eps3,hnorm,smlnum,ulp,unfl,wki,wkr
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors, and standardize the array select.
           m = 0
           pair = .false.
           do k = 1,n
              if (pair) then
                 pair = .false.
                 select(k) = .false.
              else
                 if (wi(k) == zero) then
                    if (select(k)) m = m + 1
                 else
                    pair = .true.
                    if (select(k) .or. select(k + 1)) then
                       select(k) = .true.
                       m = m + 2
                    end if
                 end if
              end if
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -13
           else if (mm < m) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('XHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_xlamch('SAFE MINIMUM')
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ldwork = n + 1
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ksr = 1
           loop_120: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are zero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == zero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == zero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_xlanhs('I',kr - kl + 1,h(kl,kl),ldh,work)
                    if (la_xisnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > zero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wkr = wr(k)
                 wki = wi(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. abs(wr(i) - wkr) + abs(wi(i) - wki) < eps3) &
                              then
                       wkr = wkr + eps3
                       go to 60
                    end if
                 end do
                 wr(k) = wkr
                 pair = wki /= zero
                 if (pair) then
                    ksi = ksr + 1
                 else
                    ksi = ksr
                 end if
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_xlaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wkr,wki,vl( &
                    kl,ksr),vl(kl,ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum, &
                              iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifaill(ksr) = k
                       ifaill(ksi) = k
                    else
                       ifaill(ksr) = 0
                       ifaill(ksi) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = 1,kl - 1
                          vl(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_xlaein(.true.,noinit,kr,h,ldh,wkr,wki,vr(1,ksr),vr(1, &
                              ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum,iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifailr(ksr) = k
                       ifailr(ksi) = k
                    else
                       ifailr(ksr) = 0
                       ifailr(ksi) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = kr + 1,n
                          vr(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (pair) then
                    ksr = ksr + 2
                 else
                    ksr = ksr + 1
                 end if
              end if
           end do loop_120
           return
     end subroutine la_xhsein
#endif
#ifdef LA_WITH_QP
     !> QHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a real upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_qhsein(side,eigsrc,initv,select,n,h,ldh,wr,wi,vl,ldvl,vr,ldvr, &
               mm,m,work,ifaill,ifailr,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(inout) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(qp),intent(in) :: h(ldh,*),wi(*)
           real(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),wr(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,pair,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ksi,ksr,ldwork
           real(qp) :: bignum,eps3,hnorm,smlnum,ulp,unfl,wki,wkr
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors, and standardize the array select.
           m = 0
           pair = .false.
           do k = 1,n
              if (pair) then
                 pair = .false.
                 select(k) = .false.
              else
                 if (wi(k) == zero) then
                    if (select(k)) m = m + 1
                 else
                    pair = .true.
                    if (select(k) .or. select(k + 1)) then
                       select(k) = .true.
                       m = m + 2
                    end if
                 end if
              end if
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -11
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -13
           else if (mm < m) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('QHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_qlamch('SAFE MINIMUM')
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           bignum = (one - ulp)/smlnum
           ldwork = n + 1
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ksr = 1
           loop_120: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are zero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == zero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == zero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_qlanhs('I',kr - kl + 1,h(kl,kl),ldh,work)
                    if (la_qisnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > zero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wkr = wr(k)
                 wki = wi(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. abs(wr(i) - wkr) + abs(wi(i) - wki) < eps3) &
                              then
                       wkr = wkr + eps3
                       go to 60
                    end if
                 end do
                 wr(k) = wkr
                 pair = wki /= zero
                 if (pair) then
                    ksi = ksr + 1
                 else
                    ksi = ksr
                 end if
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_qlaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wkr,wki,vl( &
                    kl,ksr),vl(kl,ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum, &
                              iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifaill(ksr) = k
                       ifaill(ksi) = k
                    else
                       ifaill(ksr) = 0
                       ifaill(ksi) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = 1,kl - 1
                          vl(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_qlaein(.true.,noinit,kr,h,ldh,wkr,wki,vr(1,ksr),vr(1, &
                              ksi),work,ldwork,work(n*n + n + 1),eps3,smlnum,bignum,iinfo)
                    if (iinfo > 0) then
                       if (pair) then
                          info = info + 2
                       else
                          info = info + 1
                       end if
                       ifailr(ksr) = k
                       ifailr(ksi) = k
                    else
                       ifailr(ksr) = 0
                       ifailr(ksi) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ksr) = zero
                    end do
                    if (pair) then
                       do i = kr + 1,n
                          vr(i,ksi) = zero
                       end do
                    end if
                 end if
                 if (pair) then
                    ksr = ksr + 2
                 else
                    ksr = ksr + 1
                 end if
              end if
           end do loop_120
           return
     end subroutine la_qhsein
#endif

     !> STRSEN: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
     !> the leading diagonal blocks of the upper quasi-triangular matrix T,
     !> and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.
     !> T must be in Schur canonical form (as returned by SHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_strsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work, &
               lwork,iwork,liwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,liwork,lwork,n
           real(sp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(sp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,kk,ks,liwmin,lwmin,n1,n2,nn
           real(sp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else
              ! set m to the dimension of the specified invariant subspace,
              ! and test lwork and liwork.
              m = 0
              pair = .false.
              do k = 1,n
                 if (pair) then
                    pair = .false.
                 else
                    if (k < n) then
                       if (t(k + 1,k) == zero) then
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
              n1 = m
              n2 = n - m
              nn = n1*n2
              if (wantsp) then
                 lwmin = max(1,2*nn)
                 liwmin = max(1,nn)
              else if (la_lsame(job,'N')) then
                 lwmin = max(1,n)
                 liwmin = 1
              else if (la_lsame(job,'E')) then
                 lwmin = max(1,nn)
                 liwmin = 1
              end if
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -17
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
           end if
           if (info /= 0) then
              call la_xerbla('STRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_slange('1',n,n,t,ldt,work)
              go to 40
           end if
           ! collect the selected blocks at the top-left corner of t.
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (t(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ierr = 0
                    kk = k
                    if (k /= ks) call la_strexc(compq,n,t,ldt,q,ldq,kk,ks,work,ierr)

                    if (ierr == 1 .or. ierr == 2) then
                       ! blocks too close to swap: exit.
                       info = 1
                       if (wants) s = zero
                       if (wantsp) sep = zero
                       go to 40
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_20
           if (wants) then
              ! solve sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_slacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_strsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_slange('F',n1,n2,work,n1,work)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_slacn2(nn,work(nn + 1),work,iwork,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve  t11*r - r*t22 = scale*x.
                    call la_strsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**t*r - r*t22**t = scale*x.
                    call la_strsyl('T','T',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! store the output eigenvalues in wr and wi.
           do k = 1,n
              wr(k) = t(k,k)
              wi(k) = zero
           end do
           do k = 1,n - 1
              if (t(k + 1,k) /= zero) then
                 wi(k) = sqrt(abs(t(k,k + 1)))*sqrt(abs(t(k + 1,k)))
                 wi(k + 1) = -wi(k)
              end if
           end do
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_strsen
     !> DTRSEN: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
     !> the leading diagonal blocks of the upper quasi-triangular matrix T,
     !> and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.
     !> T must be in Schur canonical form (as returned by DHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_dtrsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work, &
               lwork,iwork,liwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,liwork,lwork,n
           real(dp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(dp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,kk,ks,liwmin,lwmin,n1,n2,nn
           real(dp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else
              ! set m to the dimension of the specified invariant subspace,
              ! and test lwork and liwork.
              m = 0
              pair = .false.
              do k = 1,n
                 if (pair) then
                    pair = .false.
                 else
                    if (k < n) then
                       if (t(k + 1,k) == zero) then
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
              n1 = m
              n2 = n - m
              nn = n1*n2
              if (wantsp) then
                 lwmin = max(1,2*nn)
                 liwmin = max(1,nn)
              else if (la_lsame(job,'N')) then
                 lwmin = max(1,n)
                 liwmin = 1
              else if (la_lsame(job,'E')) then
                 lwmin = max(1,nn)
                 liwmin = 1
              end if
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -17
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
           end if
           if (info /= 0) then
              call la_xerbla('DTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_dlange('1',n,n,t,ldt,work)
              go to 40
           end if
           ! collect the selected blocks at the top-left corner of t.
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (t(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ierr = 0
                    kk = k
                    if (k /= ks) call la_dtrexc(compq,n,t,ldt,q,ldq,kk,ks,work,ierr)

                    if (ierr == 1 .or. ierr == 2) then
                       ! blocks too close to swap: exit.
                       info = 1
                       if (wants) s = zero
                       if (wantsp) sep = zero
                       go to 40
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_20
           if (wants) then
              ! solve sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_dlacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_dtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_dlange('F',n1,n2,work,n1,work)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_dlacn2(nn,work(nn + 1),work,iwork,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve  t11*r - r*t22 = scale*x.
                    call la_dtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**t*r - r*t22**t = scale*x.
                    call la_dtrsyl('T','T',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! store the output eigenvalues in wr and wi.
           do k = 1,n
              wr(k) = t(k,k)
              wi(k) = zero
           end do
           do k = 1,n - 1
              if (t(k + 1,k) /= zero) then
                 wi(k) = sqrt(abs(t(k,k + 1)))*sqrt(abs(t(k + 1,k)))
                 wi(k + 1) = -wi(k)
              end if
           end do
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dtrsen
#ifdef LA_WITH_XDP
     !> XTRSEN: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
     !> the leading diagonal blocks of the upper quasi-triangular matrix T,
     !> and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.
     !> T must be in Schur canonical form (as returned by XHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_xtrsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work, &
               lwork,iwork,liwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,liwork,lwork,n
           real(xdp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(xdp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,kk,ks,liwmin,lwmin,n1,n2,nn
           real(xdp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else
              ! set m to the dimension of the specified invariant subspace,
              ! and test lwork and liwork.
              m = 0
              pair = .false.
              do k = 1,n
                 if (pair) then
                    pair = .false.
                 else
                    if (k < n) then
                       if (t(k + 1,k) == zero) then
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
              n1 = m
              n2 = n - m
              nn = n1*n2
              if (wantsp) then
                 lwmin = max(1,2*nn)
                 liwmin = max(1,nn)
              else if (la_lsame(job,'N')) then
                 lwmin = max(1,n)
                 liwmin = 1
              else if (la_lsame(job,'E')) then
                 lwmin = max(1,nn)
                 liwmin = 1
              end if
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -17
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
           end if
           if (info /= 0) then
              call la_xerbla('XTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_xlange('1',n,n,t,ldt,work)
              go to 40
           end if
           ! collect the selected blocks at the top-left corner of t.
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (t(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ierr = 0
                    kk = k
                    if (k /= ks) call la_xtrexc(compq,n,t,ldt,q,ldq,kk,ks,work,ierr)

                    if (ierr == 1 .or. ierr == 2) then
                       ! blocks too close to swap: exit.
                       info = 1
                       if (wants) s = zero
                       if (wantsp) sep = zero
                       go to 40
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_20
           if (wants) then
              ! solve sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_xlacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_xtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_xlange('F',n1,n2,work,n1,work)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_xlacn2(nn,work(nn + 1),work,iwork,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve  t11*r - r*t22 = scale*x.
                    call la_xtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**t*r - r*t22**t = scale*x.
                    call la_xtrsyl('T','T',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! store the output eigenvalues in wr and wi.
           do k = 1,n
              wr(k) = t(k,k)
              wi(k) = zero
           end do
           do k = 1,n - 1
              if (t(k + 1,k) /= zero) then
                 wi(k) = sqrt(abs(t(k,k + 1)))*sqrt(abs(t(k + 1,k)))
                 wi(k + 1) = -wi(k)
              end if
           end do
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_xtrsen
#endif
#ifdef LA_WITH_QP
     !> QTRSEN: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
     !> the leading diagonal blocks of the upper quasi-triangular matrix T,
     !> and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.
     !> T must be in Schur canonical form (as returned by QHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_qtrsen(job,compq,select,n,t,ldt,q,ldq,wr,wi,m,s,sep,work, &
               lwork,iwork,liwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,liwork,lwork,n
           real(qp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(qp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,pair,swap,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,kk,ks,liwmin,lwmin,n1,n2,nn
           real(qp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else
              ! set m to the dimension of the specified invariant subspace,
              ! and test lwork and liwork.
              m = 0
              pair = .false.
              do k = 1,n
                 if (pair) then
                    pair = .false.
                 else
                    if (k < n) then
                       if (t(k + 1,k) == zero) then
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
              n1 = m
              n2 = n - m
              nn = n1*n2
              if (wantsp) then
                 lwmin = max(1,2*nn)
                 liwmin = max(1,nn)
              else if (la_lsame(job,'N')) then
                 lwmin = max(1,n)
                 liwmin = 1
              else if (la_lsame(job,'E')) then
                 lwmin = max(1,nn)
                 liwmin = 1
              end if
              if (lwork < lwmin .and. .not. lquery) then
                 info = -15
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -17
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
           end if
           if (info /= 0) then
              call la_xerbla('QTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_qlange('1',n,n,t,ldt,work)
              go to 40
           end if
           ! collect the selected blocks at the top-left corner of t.
           ks = 0
           pair = .false.
           loop_20: do k = 1,n
              if (pair) then
                 pair = .false.
              else
                 swap = select(k)
                 if (k < n) then
                    if (t(k + 1,k) /= zero) then
                       pair = .true.
                       swap = swap .or. select(k + 1)
                    end if
                 end if
                 if (swap) then
                    ks = ks + 1
                    ! swap the k-th block to position ks.
                    ierr = 0
                    kk = k
                    if (k /= ks) call la_qtrexc(compq,n,t,ldt,q,ldq,kk,ks,work,ierr)

                    if (ierr == 1 .or. ierr == 2) then
                       ! blocks too close to swap: exit.
                       info = 1
                       if (wants) s = zero
                       if (wantsp) sep = zero
                       go to 40
                    end if
                    if (pair) ks = ks + 1
                 end if
              end if
           end do loop_20
           if (wants) then
              ! solve sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_qlacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_qtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_qlange('F',n1,n2,work,n1,work)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_qlacn2(nn,work(nn + 1),work,iwork,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve  t11*r - r*t22 = scale*x.
                    call la_qtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**t*r - r*t22**t = scale*x.
                    call la_qtrsyl('T','T',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! store the output eigenvalues in wr and wi.
           do k = 1,n
              wr(k) = t(k,k)
              wi(k) = zero
           end do
           do k = 1,n - 1
              if (t(k + 1,k) /= zero) then
                 wi(k) = sqrt(abs(t(k,k + 1)))*sqrt(abs(t(k + 1,k)))
                 wi(k + 1) = -wi(k)
              end if
           end do
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qtrsen
#endif

     !> STRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a real upper
     !> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
     !> orthogonal).
     !> T must be in Schur canonical form (as returned by SHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_strsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm,m, &
               work,ldwork,iwork,info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(out) :: s(*),sep(*),work(ldwork,*)
           real(sp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: pair,somcon,wantbh,wants,wantsp
           integer(ilp) :: i,ierr,ifst,ilst,j,k,kase,ks,n2,nn
           real(sp) :: bignum,cond,cs,delta,dumm,eps,est,lnrm,mu,prod,prod1,prod2, &
                     rnrm,scale,smlnum,sn
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(sp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
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
                          if (t(k + 1,k) == zero) then
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
              if (mm < m) then
                 info = -13
              else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ks = 0
           pair = .false.
           loop_60: do k = 1,n
              ! determine whether t(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_60
              else
                 if (k < n) pair = t(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_60
                 else
                    if (.not. select(k)) cycle loop_60
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (.not. pair) then
                    ! real eigenvalue.
                    prod = la_sdot(n,vr(1,ks),1,vl(1,ks),1)
                    rnrm = la_snrm2(n,vr(1,ks),1)
                    lnrm = la_snrm2(n,vl(1,ks),1)
                    s(ks) = abs(prod)/(rnrm*lnrm)
                 else
                    ! complex eigenvalue.
                    prod1 = la_sdot(n,vr(1,ks),1,vl(1,ks),1)
                    prod1 = prod1 + la_sdot(n,vr(1,ks + 1),1,vl(1,ks + 1),1)
                    prod2 = la_sdot(n,vl(1,ks),1,vr(1,ks + 1),1)
                    prod2 = prod2 - la_sdot(n,vl(1,ks + 1),1,vr(1,ks),1)
                    rnrm = la_slapy2(la_snrm2(n,vr(1,ks),1),la_snrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_slapy2(la_snrm2(n,vl(1,ks),1),la_snrm2(n,vl( &
                              1,ks + 1),1))
                    cond = la_slapy2(prod1,prod2)/(rnrm*lnrm)
                    s(ks) = cond
                    s(ks + 1) = cond
                 end if
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the diagonal
                 ! block beginning at t(k,k) to the (1,1) position.
                 call la_slacpy('FULL',n,n,t,ldt,work,ldwork)
                 ifst = k
                 ilst = 1
                 call la_strexc('NO Q',n,work,ldwork,dummy,1,ifst,ilst,work(1,n + 1), &
                            ierr)
                 if (ierr == 1 .or. ierr == 2) then
                    ! could not swap because blocks not well separated
                    scale = one
                    est = bignum
                 else
                    ! reordering successful
                    if (work(2,1) == zero) then
                       ! form c = t22 - lambda*i in work(2:n,2:n).
                       do i = 2,n
                          work(i,i) = work(i,i) - work(1,1)
                       end do
                       n2 = 1
                       nn = n - 1
                    else
                       ! triangularize the 2 by 2 block by unitary
                       ! transformation u = [  cs   i*ss ]
                                          ! [ i*ss   cs  ].
                       ! such that the (1,1) position of work is complex
                       ! eigenvalue lambda with positive imaginary part. (2,2)
                       ! position of work is the complex eigenvalue lambda
                       ! with negative imaginary  part.
                       mu = sqrt(abs(work(1,2)))*sqrt(abs(work(2,1)))
                       delta = la_slapy2(mu,work(2,1))
                       cs = mu/delta
                       sn = -work(2,1)/delta
                       ! form
                       ! c**t = work(2:n,2:n) + i*[rwork(1) ..... rwork(n-1) ]
                                                ! [   mu                     ]
                                                ! [         ..               ]
                                                ! [             ..           ]
                                                ! [                  mu      ]
                       ! where c**t is transpose of matrix c,
                       ! and rwork is stored starting in the n+1-st column of
                       ! work.
                       do j = 3,n
                          work(2,j) = cs*work(2,j)
                          work(j,j) = work(j,j) - work(1,1)
                       end do
                       work(2,2) = zero
                       work(1,n + 1) = two*mu
                       do i = 2,n - 1
                          work(i,n + 1) = sn*work(1,i + 1)
                       end do
                       n2 = 2
                       nn = 2*(n - 1)
                    end if
                    ! estimate norm(inv(c**t))
                    est = zero
                    kase = 0
                    50 continue
                    call la_slacn2(nn,work(1,n + 2),work(1,n + 4),iwork,est,kase, &
                              isave)
                    if (kase /= 0) then
                       if (kase == 1) then
                          if (n2 == 1) then
                             ! real eigenvalue: solve c**t*x = scale*c.
                             call la_slaqtr(.true.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                       dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c**t*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_slaqtr(.true.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       else
                          if (n2 == 1) then
                             ! real eigenvalue: solve c*x = scale*c.
                             call la_slaqtr(.false.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                        dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_slaqtr(.false.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       end if
                       go to 50
                    end if
                 end if
                 sep(ks) = scale/max(est,smlnum)
                 if (pair) sep(ks + 1) = sep(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_60
           return
     end subroutine la_strsna
     !> DTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a real upper
     !> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
     !> orthogonal).
     !> T must be in Schur canonical form (as returned by DHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_dtrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm,m, &
               work,ldwork,iwork,info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(out) :: s(*),sep(*),work(ldwork,*)
           real(dp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: pair,somcon,wantbh,wants,wantsp
           integer(ilp) :: i,ierr,ifst,ilst,j,k,kase,ks,n2,nn
           real(dp) :: bignum,cond,cs,delta,dumm,eps,est,lnrm,mu,prod,prod1,prod2, &
                     rnrm,scale,smlnum,sn
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(dp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
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
                          if (t(k + 1,k) == zero) then
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
              if (mm < m) then
                 info = -13
              else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ks = 0
           pair = .false.
           loop_60: do k = 1,n
              ! determine whether t(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_60
              else
                 if (k < n) pair = t(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_60
                 else
                    if (.not. select(k)) cycle loop_60
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (.not. pair) then
                    ! real eigenvalue.
                    prod = la_ddot(n,vr(1,ks),1,vl(1,ks),1)
                    rnrm = la_dnrm2(n,vr(1,ks),1)
                    lnrm = la_dnrm2(n,vl(1,ks),1)
                    s(ks) = abs(prod)/(rnrm*lnrm)
                 else
                    ! complex eigenvalue.
                    prod1 = la_ddot(n,vr(1,ks),1,vl(1,ks),1)
                    prod1 = prod1 + la_ddot(n,vr(1,ks + 1),1,vl(1,ks + 1),1)
                    prod2 = la_ddot(n,vl(1,ks),1,vr(1,ks + 1),1)
                    prod2 = prod2 - la_ddot(n,vl(1,ks + 1),1,vr(1,ks),1)
                    rnrm = la_dlapy2(la_dnrm2(n,vr(1,ks),1),la_dnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_dlapy2(la_dnrm2(n,vl(1,ks),1),la_dnrm2(n,vl( &
                              1,ks + 1),1))
                    cond = la_dlapy2(prod1,prod2)/(rnrm*lnrm)
                    s(ks) = cond
                    s(ks + 1) = cond
                 end if
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the diagonal
                 ! block beginning at t(k,k) to the (1,1) position.
                 call la_dlacpy('FULL',n,n,t,ldt,work,ldwork)
                 ifst = k
                 ilst = 1
                 call la_dtrexc('NO Q',n,work,ldwork,dummy,1,ifst,ilst,work(1,n + 1), &
                            ierr)
                 if (ierr == 1 .or. ierr == 2) then
                    ! could not swap because blocks not well separated
                    scale = one
                    est = bignum
                 else
                    ! reordering successful
                    if (work(2,1) == zero) then
                       ! form c = t22 - lambda*i in work(2:n,2:n).
                       do i = 2,n
                          work(i,i) = work(i,i) - work(1,1)
                       end do
                       n2 = 1
                       nn = n - 1
                    else
                       ! triangularize the 2 by 2 block by unitary
                       ! transformation u = [  cs   i*ss ]
                                          ! [ i*ss   cs  ].
                       ! such that the (1,1) position of work is complex
                       ! eigenvalue lambda with positive imaginary part. (2,2)
                       ! position of work is the complex eigenvalue lambda
                       ! with negative imaginary  part.
                       mu = sqrt(abs(work(1,2)))*sqrt(abs(work(2,1)))
                       delta = la_dlapy2(mu,work(2,1))
                       cs = mu/delta
                       sn = -work(2,1)/delta
                       ! form
                       ! c**t = work(2:n,2:n) + i*[rwork(1) ..... rwork(n-1) ]
                                                ! [   mu                     ]
                                                ! [         ..               ]
                                                ! [             ..           ]
                                                ! [                  mu      ]
                       ! where c**t is transpose of matrix c,
                       ! and rwork is stored starting in the n+1-st column of
                       ! work.
                       do j = 3,n
                          work(2,j) = cs*work(2,j)
                          work(j,j) = work(j,j) - work(1,1)
                       end do
                       work(2,2) = zero
                       work(1,n + 1) = two*mu
                       do i = 2,n - 1
                          work(i,n + 1) = sn*work(1,i + 1)
                       end do
                       n2 = 2
                       nn = 2*(n - 1)
                    end if
                    ! estimate norm(inv(c**t))
                    est = zero
                    kase = 0
                    50 continue
                    call la_dlacn2(nn,work(1,n + 2),work(1,n + 4),iwork,est,kase, &
                              isave)
                    if (kase /= 0) then
                       if (kase == 1) then
                          if (n2 == 1) then
                             ! real eigenvalue: solve c**t*x = scale*c.
                             call la_dlaqtr(.true.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                       dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c**t*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_dlaqtr(.true.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       else
                          if (n2 == 1) then
                             ! real eigenvalue: solve c*x = scale*c.
                             call la_dlaqtr(.false.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                        dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_dlaqtr(.false.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       end if
                       go to 50
                    end if
                 end if
                 sep(ks) = scale/max(est,smlnum)
                 if (pair) sep(ks + 1) = sep(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_60
           return
     end subroutine la_dtrsna
#ifdef LA_WITH_XDP
     !> XTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a real upper
     !> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
     !> orthogonal).
     !> T must be in Schur canonical form (as returned by XHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_xtrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm,m, &
               work,ldwork,iwork,info)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(out) :: s(*),sep(*),work(ldwork,*)
           real(xdp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: pair,somcon,wantbh,wants,wantsp
           integer(ilp) :: i,ierr,ifst,ilst,j,k,kase,ks,n2,nn
           real(xdp) :: bignum,cond,cs,delta,dumm,eps,est,lnrm,mu,prod,prod1,prod2, &
                     rnrm,scale,smlnum,sn
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(xdp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
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
                          if (t(k + 1,k) == zero) then
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
              if (mm < m) then
                 info = -13
              else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ks = 0
           pair = .false.
           loop_60: do k = 1,n
              ! determine whether t(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_60
              else
                 if (k < n) pair = t(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_60
                 else
                    if (.not. select(k)) cycle loop_60
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (.not. pair) then
                    ! real eigenvalue.
                    prod = la_xdot(n,vr(1,ks),1,vl(1,ks),1)
                    rnrm = la_xnrm2(n,vr(1,ks),1)
                    lnrm = la_xnrm2(n,vl(1,ks),1)
                    s(ks) = abs(prod)/(rnrm*lnrm)
                 else
                    ! complex eigenvalue.
                    prod1 = la_xdot(n,vr(1,ks),1,vl(1,ks),1)
                    prod1 = prod1 + la_xdot(n,vr(1,ks + 1),1,vl(1,ks + 1),1)
                    prod2 = la_xdot(n,vl(1,ks),1,vr(1,ks + 1),1)
                    prod2 = prod2 - la_xdot(n,vl(1,ks + 1),1,vr(1,ks),1)
                    rnrm = la_xlapy2(la_xnrm2(n,vr(1,ks),1),la_xnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_xlapy2(la_xnrm2(n,vl(1,ks),1),la_xnrm2(n,vl( &
                              1,ks + 1),1))
                    cond = la_xlapy2(prod1,prod2)/(rnrm*lnrm)
                    s(ks) = cond
                    s(ks + 1) = cond
                 end if
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the diagonal
                 ! block beginning at t(k,k) to the (1,1) position.
                 call la_xlacpy('FULL',n,n,t,ldt,work,ldwork)
                 ifst = k
                 ilst = 1
                 call la_xtrexc('NO Q',n,work,ldwork,dummy,1,ifst,ilst,work(1,n + 1), &
                            ierr)
                 if (ierr == 1 .or. ierr == 2) then
                    ! could not swap because blocks not well separated
                    scale = one
                    est = bignum
                 else
                    ! reordering successful
                    if (work(2,1) == zero) then
                       ! form c = t22 - lambda*i in work(2:n,2:n).
                       do i = 2,n
                          work(i,i) = work(i,i) - work(1,1)
                       end do
                       n2 = 1
                       nn = n - 1
                    else
                       ! triangularize the 2 by 2 block by unitary
                       ! transformation u = [  cs   i*ss ]
                                          ! [ i*ss   cs  ].
                       ! such that the (1,1) position of work is complex
                       ! eigenvalue lambda with positive imaginary part. (2,2)
                       ! position of work is the complex eigenvalue lambda
                       ! with negative imaginary  part.
                       mu = sqrt(abs(work(1,2)))*sqrt(abs(work(2,1)))
                       delta = la_xlapy2(mu,work(2,1))
                       cs = mu/delta
                       sn = -work(2,1)/delta
                       ! form
                       ! c**t = work(2:n,2:n) + i*[rwork(1) ..... rwork(n-1) ]
                                                ! [   mu                     ]
                                                ! [         ..               ]
                                                ! [             ..           ]
                                                ! [                  mu      ]
                       ! where c**t is transpose of matrix c,
                       ! and rwork is stored starting in the n+1-st column of
                       ! work.
                       do j = 3,n
                          work(2,j) = cs*work(2,j)
                          work(j,j) = work(j,j) - work(1,1)
                       end do
                       work(2,2) = zero
                       work(1,n + 1) = two*mu
                       do i = 2,n - 1
                          work(i,n + 1) = sn*work(1,i + 1)
                       end do
                       n2 = 2
                       nn = 2*(n - 1)
                    end if
                    ! estimate norm(inv(c**t))
                    est = zero
                    kase = 0
                    50 continue
                    call la_xlacn2(nn,work(1,n + 2),work(1,n + 4),iwork,est,kase, &
                              isave)
                    if (kase /= 0) then
                       if (kase == 1) then
                          if (n2 == 1) then
                             ! real eigenvalue: solve c**t*x = scale*c.
                             call la_xlaqtr(.true.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                       dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c**t*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_xlaqtr(.true.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       else
                          if (n2 == 1) then
                             ! real eigenvalue: solve c*x = scale*c.
                             call la_xlaqtr(.false.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                        dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_xlaqtr(.false.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       end if
                       go to 50
                    end if
                 end if
                 sep(ks) = scale/max(est,smlnum)
                 if (pair) sep(ks + 1) = sep(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_60
           return
     end subroutine la_xtrsna
#endif
#ifdef LA_WITH_QP
     !> QTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a real upper
     !> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
     !> orthogonal).
     !> T must be in Schur canonical form (as returned by QHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_qtrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm,m, &
               work,ldwork,iwork,info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(out) :: s(*),sep(*),work(ldwork,*)
           real(qp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: pair,somcon,wantbh,wants,wantsp
           integer(ilp) :: i,ierr,ifst,ilst,j,k,kase,ks,n2,nn
           real(qp) :: bignum,cond,cs,delta,dumm,eps,est,lnrm,mu,prod,prod1,prod2, &
                     rnrm,scale,smlnum,sn
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(qp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
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
                          if (t(k + 1,k) == zero) then
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
              if (mm < m) then
                 info = -13
              else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
                 info = -16
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ks = 0
           pair = .false.
           loop_60: do k = 1,n
              ! determine whether t(k,k) begins a 1-by-1 or 2-by-2 block.
              if (pair) then
                 pair = .false.
                 cycle loop_60
              else
                 if (k < n) pair = t(k + 1,k) /= zero
              end if
              ! determine whether condition numbers are required for the k-th
              ! eigenpair.
              if (somcon) then
                 if (pair) then
                    if (.not. select(k) .and. .not. select(k + 1)) cycle loop_60
                 else
                    if (.not. select(k)) cycle loop_60
                 end if
              end if
              ks = ks + 1
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 if (.not. pair) then
                    ! real eigenvalue.
                    prod = la_qdot(n,vr(1,ks),1,vl(1,ks),1)
                    rnrm = la_qnrm2(n,vr(1,ks),1)
                    lnrm = la_qnrm2(n,vl(1,ks),1)
                    s(ks) = abs(prod)/(rnrm*lnrm)
                 else
                    ! complex eigenvalue.
                    prod1 = la_qdot(n,vr(1,ks),1,vl(1,ks),1)
                    prod1 = prod1 + la_qdot(n,vr(1,ks + 1),1,vl(1,ks + 1),1)
                    prod2 = la_qdot(n,vl(1,ks),1,vr(1,ks + 1),1)
                    prod2 = prod2 - la_qdot(n,vl(1,ks + 1),1,vr(1,ks),1)
                    rnrm = la_qlapy2(la_qnrm2(n,vr(1,ks),1),la_qnrm2(n,vr( &
                              1,ks + 1),1))
                    lnrm = la_qlapy2(la_qnrm2(n,vl(1,ks),1),la_qnrm2(n,vl( &
                              1,ks + 1),1))
                    cond = la_qlapy2(prod1,prod2)/(rnrm*lnrm)
                    s(ks) = cond
                    s(ks + 1) = cond
                 end if
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the diagonal
                 ! block beginning at t(k,k) to the (1,1) position.
                 call la_qlacpy('FULL',n,n,t,ldt,work,ldwork)
                 ifst = k
                 ilst = 1
                 call la_qtrexc('NO Q',n,work,ldwork,dummy,1,ifst,ilst,work(1,n + 1), &
                            ierr)
                 if (ierr == 1 .or. ierr == 2) then
                    ! could not swap because blocks not well separated
                    scale = one
                    est = bignum
                 else
                    ! reordering successful
                    if (work(2,1) == zero) then
                       ! form c = t22 - lambda*i in work(2:n,2:n).
                       do i = 2,n
                          work(i,i) = work(i,i) - work(1,1)
                       end do
                       n2 = 1
                       nn = n - 1
                    else
                       ! triangularize the 2 by 2 block by unitary
                       ! transformation u = [  cs   i*ss ]
                                          ! [ i*ss   cs  ].
                       ! such that the (1,1) position of work is complex
                       ! eigenvalue lambda with positive imaginary part. (2,2)
                       ! position of work is the complex eigenvalue lambda
                       ! with negative imaginary  part.
                       mu = sqrt(abs(work(1,2)))*sqrt(abs(work(2,1)))
                       delta = la_qlapy2(mu,work(2,1))
                       cs = mu/delta
                       sn = -work(2,1)/delta
                       ! form
                       ! c**t = work(2:n,2:n) + i*[rwork(1) ..... rwork(n-1) ]
                                                ! [   mu                     ]
                                                ! [         ..               ]
                                                ! [             ..           ]
                                                ! [                  mu      ]
                       ! where c**t is transpose of matrix c,
                       ! and rwork is stored starting in the n+1-st column of
                       ! work.
                       do j = 3,n
                          work(2,j) = cs*work(2,j)
                          work(j,j) = work(j,j) - work(1,1)
                       end do
                       work(2,2) = zero
                       work(1,n + 1) = two*mu
                       do i = 2,n - 1
                          work(i,n + 1) = sn*work(1,i + 1)
                       end do
                       n2 = 2
                       nn = 2*(n - 1)
                    end if
                    ! estimate norm(inv(c**t))
                    est = zero
                    kase = 0
                    50 continue
                    call la_qlacn2(nn,work(1,n + 2),work(1,n + 4),iwork,est,kase, &
                              isave)
                    if (kase /= 0) then
                       if (kase == 1) then
                          if (n2 == 1) then
                             ! real eigenvalue: solve c**t*x = scale*c.
                             call la_qlaqtr(.true.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                       dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c**t*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_qlaqtr(.true.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       else
                          if (n2 == 1) then
                             ! real eigenvalue: solve c*x = scale*c.
                             call la_qlaqtr(.false.,.true.,n - 1,work(2,2),ldwork,dummy, &
                                        dumm,scale,work(1,n + 4),work(1,n + 6),ierr)
                          else
                             ! complex eigenvalue: solve
                             ! c*(p+iq) = scale*(c+id) in real arithmetic.
                             call la_qlaqtr(.false.,.false.,n - 1,work(2,2),ldwork,work( &
                                       1,n + 1),mu,scale,work(1,n + 4),work(1,n + 6),ierr)
                          end if
                       end if
                       go to 50
                    end if
                 end if
                 sep(ks) = scale/max(est,smlnum)
                 if (pair) sep(ks + 1) = sep(ks)
              end if
              if (pair) ks = ks + 1
           end do loop_60
           return
     end subroutine la_qtrsna
#endif

     !> SHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**T, where T is an upper quasi-triangular matrix (the
     !> Schur form), and Z is the orthogonal matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input orthogonal
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.

     subroutine la_shseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           real(sp),intent(inout) :: h(ldh,*),z(ldz,*)
           real(sp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_slahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_slahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           real(sp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: i,kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: max,min,real
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = real(max(1,n),KIND=sp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -11
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -13
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('SHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_slaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                        work,lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=sp),work(1))
              return
           else
              ! ==== copy eigenvalues isolated by la_sgebal ====
              do i = 1,ilo - 1
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              do i = ihi + 1,n
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              ! ==== initialize z, if requested ====
              if (initz) call la_slaset('A',n,n,zero,one,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 wr(ilo) = h(ilo,ilo)
                 wi(ilo) = zero
                 return
              end if
              ! ==== la_slahqr/la_slaqr0 crossover point ====
              nmin = la_ilaenv(12,'SHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_slaqr0 for big matrices; la_slahqr for small ones ====
              if (n > nmin) then
                 call la_slaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           work,lwork,info)
              else
                 ! ==== small matrix ====
                 call la_slahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           info)
                 if (info > 0) then
                    ! ==== a rare la_slahqr failure!  la_slaqr0 sometimes succeeds
                    ! .    when la_slahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_slaqr0 directly. ====
                       call la_slaqr0(wantt,wantz,n,ilo,kbot,h,ldh,wr,wi,ilo,ihi,z, &
                                  ldz,work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_slaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_slaqr0. ====
                       call la_slacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = zero
                       call la_slaset('A',nl,nl - n,zero,zero,hl(1,n + 1),nl)
                       call la_slaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,wr,wi,ilo,ihi, &
                                 z,ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_slacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_slaset('L',n - 2,n - 2,zero,zero, &
                         h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=sp),work(1))
           end if
     end subroutine la_shseqr
     !> DHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**T, where T is an upper quasi-triangular matrix (the
     !> Schur form), and Z is the orthogonal matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input orthogonal
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.

     subroutine la_dhseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           real(dp),intent(inout) :: h(ldh,*),z(ldz,*)
           real(dp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_dlahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_dlahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           real(dp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: i,kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = real(max(1,n),KIND=dp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -11
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -13
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('DHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_dlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                        work,lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=dp),work(1))
              return
           else
              ! ==== copy eigenvalues isolated by la_dgebal ====
              do i = 1,ilo - 1
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              do i = ihi + 1,n
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              ! ==== initialize z, if requested ====
              if (initz) call la_dlaset('A',n,n,zero,one,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 wr(ilo) = h(ilo,ilo)
                 wi(ilo) = zero
                 return
              end if
              ! ==== la_dlahqr/la_dlaqr0 crossover point ====
              nmin = la_ilaenv(12,'DHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_dlaqr0 for big matrices; la_dlahqr for small ones ====
              if (n > nmin) then
                 call la_dlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           work,lwork,info)
              else
                 ! ==== small matrix ====
                 call la_dlahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           info)
                 if (info > 0) then
                    ! ==== a rare la_dlahqr failure!  la_dlaqr0 sometimes succeeds
                    ! .    when la_dlahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_dlaqr0 directly. ====
                       call la_dlaqr0(wantt,wantz,n,ilo,kbot,h,ldh,wr,wi,ilo,ihi,z, &
                                  ldz,work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_dlaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_dlaqr0. ====
                       call la_dlacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = zero
                       call la_dlaset('A',nl,nl - n,zero,zero,hl(1,n + 1),nl)
                       call la_dlaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,wr,wi,ilo,ihi, &
                                 z,ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_dlacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_dlaset('L',n - 2,n - 2,zero,zero, &
                         h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=dp),work(1))
           end if
     end subroutine la_dhseqr
#ifdef LA_WITH_XDP
     !> XHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**T, where T is an upper quasi-triangular matrix (the
     !> Schur form), and Z is the orthogonal matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input orthogonal
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.

     subroutine la_xhseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           real(xdp),intent(inout) :: h(ldh,*),z(ldz,*)
           real(xdp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_xlahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_xlahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           real(xdp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: i,kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = real(max(1,n),KIND=xdp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -11
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -13
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('XHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_xlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                        work,lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=xdp),work(1))
              return
           else
              ! ==== copy eigenvalues isolated by la_xgebal ====
              do i = 1,ilo - 1
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              do i = ihi + 1,n
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              ! ==== initialize z, if requested ====
              if (initz) call la_xlaset('A',n,n,zero,one,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 wr(ilo) = h(ilo,ilo)
                 wi(ilo) = zero
                 return
              end if
              ! ==== la_xlahqr/la_xlaqr0 crossover point ====
              nmin = la_ilaenv(12,'XHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_xlaqr0 for big matrices; la_xlahqr for small ones ====
              if (n > nmin) then
                 call la_xlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           work,lwork,info)
              else
                 ! ==== small matrix ====
                 call la_xlahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           info)
                 if (info > 0) then
                    ! ==== a rare la_xlahqr failure!  la_xlaqr0 sometimes succeeds
                    ! .    when la_xlahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_xlaqr0 directly. ====
                       call la_xlaqr0(wantt,wantz,n,ilo,kbot,h,ldh,wr,wi,ilo,ihi,z, &
                                  ldz,work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_xlaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_xlaqr0. ====
                       call la_xlacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = zero
                       call la_xlaset('A',nl,nl - n,zero,zero,hl(1,n + 1),nl)
                       call la_xlaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,wr,wi,ilo,ihi, &
                                 z,ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_xlacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_xlaset('L',n - 2,n - 2,zero,zero, &
                         h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=xdp),work(1))
           end if
     end subroutine la_xhseqr
#endif
#ifdef LA_WITH_QP
     !> QHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**T, where T is an upper quasi-triangular matrix (the
     !> Schur form), and Z is the orthogonal matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input orthogonal
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.

     subroutine la_qhseqr(job,compz,n,ilo,ihi,h,ldh,wr,wi,z,ldz,work,lwork,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           real(qp),intent(inout) :: h(ldh,*),z(ldz,*)
           real(qp),intent(out) :: wi(*),work(*),wr(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_qlahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_qlahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           real(qp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: i,kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = real(max(1,n),KIND=qp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -11
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -13
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('QHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_qlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                        work,lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=qp),work(1))
              return
           else
              ! ==== copy eigenvalues isolated by la_qgebal ====
              do i = 1,ilo - 1
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              do i = ihi + 1,n
                 wr(i) = h(i,i)
                 wi(i) = zero
              end do
              ! ==== initialize z, if requested ====
              if (initz) call la_qlaset('A',n,n,zero,one,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 wr(ilo) = h(ilo,ilo)
                 wi(ilo) = zero
                 return
              end if
              ! ==== la_qlahqr/la_qlaqr0 crossover point ====
              nmin = la_ilaenv(12,'QHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_qlaqr0 for big matrices; la_qlahqr for small ones ====
              if (n > nmin) then
                 call la_qlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           work,lwork,info)
              else
                 ! ==== small matrix ====
                 call la_qlahqr(wantt,wantz,n,ilo,ihi,h,ldh,wr,wi,ilo,ihi,z,ldz, &
                           info)
                 if (info > 0) then
                    ! ==== a rare la_qlahqr failure!  la_qlaqr0 sometimes succeeds
                    ! .    when la_qlahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_qlaqr0 directly. ====
                       call la_qlaqr0(wantt,wantz,n,ilo,kbot,h,ldh,wr,wi,ilo,ihi,z, &
                                  ldz,work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_qlaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_qlaqr0. ====
                       call la_qlacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = zero
                       call la_qlaset('A',nl,nl - n,zero,zero,hl(1,n + 1),nl)
                       call la_qlaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,wr,wi,ilo,ihi, &
                                 z,ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_qlacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_qlaset('L',n - 2,n - 2,zero,zero, &
                         h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = max(real(max(1,n),KIND=qp),work(1))
           end if
     end subroutine la_qhseqr
#endif

     !> CTREVC: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.

     pure subroutine la_ctrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           complex(sp),parameter :: cmzero = (0.0e+0_sp,0.0e+0_sp)
           complex(sp),parameter :: cmone = (1.0e+0_sp,0.0e+0_sp)

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki
           real(sp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(sp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_slamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_slabad(unfl,ovfl)
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i + n) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_scasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! compute right eigenvectors.
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(1) = cmone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k) = -t(k,ki)
                 end do
                 ! solve the triangular system:
                    ! (t(1:ki-1,1:ki-1) - t(ki,ki))*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_clatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    call la_ccopy(ki,work(1),1,vr(1,is),1)
                    ii = la_icamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_csscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = cmzero
                    end do
                 else
                    if (ki > 1) call la_cgemv('N',n,ki - 1,cmone,vr,ldvr,work(1),1, &
                              cmplx(scale,KIND=sp),vr(1,ki),1)
                    ii = la_icamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_csscal(n,remax,vr(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k + n)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! compute left eigenvectors.
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(n) = cmone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k) = -conjg(t(ki,k))
                 end do
                 ! solve the triangular system:
                    ! (t(ki+1:n,ki+1:n) - t(ki,ki))**h*x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_clatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    call la_ccopy(n - ki + 1,work(ki),1,vl(ki,is),1)
                    ii = la_icamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_csscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = cmzero
                    end do
                 else
                    if (ki < n) call la_cgemv('N',n,n - ki,cmone,vl(1,ki + 1),ldvl,work( &
                              ki + 1),1,cmplx(scale,KIND=sp),vl(1,ki),1)
                    ii = la_icamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_csscal(n,remax,vl(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k + n)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ctrevc
     !> ZTREVC: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.

     pure subroutine la_ztrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           complex(dp),parameter :: cmzero = (0.0e+0_dp,0.0e+0_dp)
           complex(dp),parameter :: cmone = (1.0e+0_dp,0.0e+0_dp)

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki
           real(dp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(dp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_dlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_dlabad(unfl,ovfl)
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i + n) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_dzasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! compute right eigenvectors.
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(1) = cmone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k) = -t(k,ki)
                 end do
                 ! solve the triangular system:
                    ! (t(1:ki-1,1:ki-1) - t(ki,ki))*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_zlatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    call la_zcopy(ki,work(1),1,vr(1,is),1)
                    ii = la_izamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_zdscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = cmzero
                    end do
                 else
                    if (ki > 1) call la_zgemv('N',n,ki - 1,cmone,vr,ldvr,work(1),1, &
                              cmplx(scale,KIND=dp),vr(1,ki),1)
                    ii = la_izamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_zdscal(n,remax,vr(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k + n)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! compute left eigenvectors.
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(n) = cmone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k) = -conjg(t(ki,k))
                 end do
                 ! solve the triangular system:
                    ! (t(ki+1:n,ki+1:n) - t(ki,ki))**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_zlatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    call la_zcopy(n - ki + 1,work(ki),1,vl(ki,is),1)
                    ii = la_izamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_zdscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = cmzero
                    end do
                 else
                    if (ki < n) call la_zgemv('N',n,n - ki,cmone,vl(1,ki + 1),ldvl,work( &
                              ki + 1),1,cmplx(scale,KIND=dp),vl(1,ki),1)
                    ii = la_izamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_zdscal(n,remax,vl(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k + n)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ztrevc
#ifdef LA_WITH_XDP
     !> YTREVC: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by YHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.

     pure subroutine la_ytrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,rwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           complex(xdp),parameter :: cmzero = (0.0e+0_xdp,0.0e+0_xdp)
           complex(xdp),parameter :: cmone = (1.0e+0_xdp,0.0e+0_xdp)

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki
           real(xdp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(xdp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('YTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_xlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_xlabad(unfl,ovfl)
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i + n) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_xyasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! compute right eigenvectors.
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(1) = cmone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k) = -t(k,ki)
                 end do
                 ! solve the triangular system:
                    ! (t(1:ki-1,1:ki-1) - t(ki,ki))*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_ylatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    call la_ycopy(ki,work(1),1,vr(1,is),1)
                    ii = la_iyamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_yxscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = cmzero
                    end do
                 else
                    if (ki > 1) call la_ygemv('N',n,ki - 1,cmone,vr,ldvr,work(1),1, &
                              cmplx(scale,KIND=xdp),vr(1,ki),1)
                    ii = la_iyamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_yxscal(n,remax,vr(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k + n)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! compute left eigenvectors.
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(n) = cmone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k) = -conjg(t(ki,k))
                 end do
                 ! solve the triangular system:
                    ! (t(ki+1:n,ki+1:n) - t(ki,ki))**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_ylatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    call la_ycopy(n - ki + 1,work(ki),1,vl(ki,is),1)
                    ii = la_iyamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_yxscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = cmzero
                    end do
                 else
                    if (ki < n) call la_ygemv('N',n,n - ki,cmone,vl(1,ki + 1),ldvl,work( &
                              ki + 1),1,cmplx(scale,KIND=xdp),vl(1,ki),1)
                    ii = la_iyamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_yxscal(n,remax,vl(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k + n)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ytrevc
#endif
#ifdef LA_WITH_QP
     !> WTREVC: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by WHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix.  If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.

     pure subroutine la_wtrevc(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           complex(qp),parameter :: cmzero = (0.0e+0_qp,0.0e+0_qp)
           complex(qp),parameter :: cmone = (1.0e+0_qp,0.0e+0_qp)

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki
           real(qp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(qp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WTREVC',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set the constants to control overflow.
           unfl = la_qlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_qlabad(unfl,ovfl)
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i + n) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_qwasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! compute right eigenvectors.
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(1) = cmone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k) = -t(k,ki)
                 end do
                 ! solve the triangular system:
                    ! (t(1:ki-1,1:ki-1) - t(ki,ki))*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_wlatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    call la_wcopy(ki,work(1),1,vr(1,is),1)
                    ii = la_iwamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_wqscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = cmzero
                    end do
                 else
                    if (ki > 1) call la_wgemv('N',n,ki - 1,cmone,vr,ldvr,work(1),1, &
                              cmplx(scale,KIND=qp),vr(1,ki),1)
                    ii = la_iwamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_wqscal(n,remax,vr(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k + n)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! compute left eigenvectors.
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 work(n) = cmone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k) = -conjg(t(ki,k))
                 end do
                 ! solve the triangular system:
                    ! (t(ki+1:n,ki+1:n) - t(ki,ki))**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_wlatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1),scale,rwork,info)
                    work(ki) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    call la_wcopy(n - ki + 1,work(ki),1,vl(ki,is),1)
                    ii = la_iwamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_wqscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = cmzero
                    end do
                 else
                    if (ki < n) call la_wgemv('N',n,n - ki,cmone,vl(1,ki + 1),ldvl,work( &
                              ki + 1),1,cmplx(scale,KIND=qp),vl(1,ki),1)
                    ii = la_iwamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_wqscal(n,remax,vl(1,ki),1)
                 end if
                 ! set back the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k + n)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_wtrevc
#endif

     !> CTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_ctrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,rwork,lrwork,info)
        use la_constants_sp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,lrwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki,iv,maxwrk,nb
           real(sp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(sp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           nb = la_ilaenv(1,'CTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           rwork(1) = n
           lquery = (lwork == -1 .or. lrwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -14
           else if (lrwork < max(1,n) .and. .not. lquery) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('CTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_claset('F',n,1 + 2*nb,czero,czero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_slamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_slabad(unfl,ovfl)
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_scasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=nb=1;
              ! blocked     version starts with iv=nb, goes down to 1.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = nb
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex right eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k + iv*n) = -t(k,ki)
                 end do
                 ! solve upper triangular system:
                 ! [ t(1:ki-1,1:ki-1) - t(ki,ki) ]*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_clatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vr and normalize.
                    call la_ccopy(ki,work(1 + iv*n),1,vr(1,is),1)
                    ii = la_icamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_csscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki > 1) call la_cgemv('N',n,ki - 1,cone,vr,ldvr,work(1 + iv*n),1, &
                               cmplx(scale,KIND=sp),vr(1,ki),1)
                    ii = la_icamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_csscal(n,remax,vr(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out below vector
                    do k = ki + 1,n
                       work(k + iv*n) = czero
                    end do
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == 1) .or. (ki == 1)) then
                       call la_cgemm('N','N',n,nb - iv + 1,ki + nb - iv,cone,vr,ldvr,work(1 + &
                                 (iv)*n),n,czero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          ii = la_icamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_csscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_clacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = 1
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex left eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k + iv*n) = -conjg(t(ki,k))
                 end do
                 ! solve conjugate-transposed triangular system:
                 ! [ t(ki+1:n,ki+1:n) - t(ki,ki) ]**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_clatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vl and normalize.
                    call la_ccopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                    ii = la_icamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_csscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki < n) call la_cgemv('N',n,n - ki,cone,vl(1,ki + 1),ldvl,work(ki + &
                              1 + iv*n),1,cmplx(scale,KIND=sp),vl(1,ki),1)
                    ii = la_icamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_csscal(n,remax,vl(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out above vector
                    ! could go from ki-nv+1 to ki-1
                    do k = 1,ki - 1
                       work(k + iv*n) = czero
                    end do
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == nb) .or. (ki == n)) then
                       call la_cgemm('N','N',n,iv,n - ki + iv,cone,vl(1,ki - iv + 1),ldvl, &
                                 work(ki - iv + 1 + (1)*n),n,czero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          ii = la_icamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_csscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_clacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ctrevc3
     !> ZTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_ztrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,rwork,lrwork,info)
        use la_constants_dp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,lrwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki,iv,maxwrk,nb
           real(dp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(dp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           nb = la_ilaenv(1,'ZTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           rwork(1) = n
           lquery = (lwork == -1 .or. lrwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -14
           else if (lrwork < max(1,n) .and. .not. lquery) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('ZTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_zlaset('F',n,1 + 2*nb,czero,czero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_dlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_dlabad(unfl,ovfl)
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_dzasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=nb=1;
              ! blocked     version starts with iv=nb, goes down to 1.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = nb
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex right eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k + iv*n) = -t(k,ki)
                 end do
                 ! solve upper triangular system:
                 ! [ t(1:ki-1,1:ki-1) - t(ki,ki) ]*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_zlatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vr and normalize.
                    call la_zcopy(ki,work(1 + iv*n),1,vr(1,is),1)
                    ii = la_izamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_zdscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki > 1) call la_zgemv('N',n,ki - 1,cone,vr,ldvr,work(1 + iv*n),1, &
                               cmplx(scale,KIND=dp),vr(1,ki),1)
                    ii = la_izamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_zdscal(n,remax,vr(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out below vector
                    do k = ki + 1,n
                       work(k + iv*n) = czero
                    end do
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == 1) .or. (ki == 1)) then
                       call la_zgemm('N','N',n,nb - iv + 1,ki + nb - iv,cone,vr,ldvr,work(1 + &
                                 (iv)*n),n,czero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          ii = la_izamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_zdscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_zlacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = 1
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex left eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k + iv*n) = -conjg(t(ki,k))
                 end do
                 ! solve conjugate-transposed triangular system:
                 ! [ t(ki+1:n,ki+1:n) - t(ki,ki) ]**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_zlatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vl and normalize.
                    call la_zcopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                    ii = la_izamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_zdscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki < n) call la_zgemv('N',n,n - ki,cone,vl(1,ki + 1),ldvl,work(ki + &
                              1 + iv*n),1,cmplx(scale,KIND=dp),vl(1,ki),1)
                    ii = la_izamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_zdscal(n,remax,vl(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out above vector
                    ! could go from ki-nv+1 to ki-1
                    do k = 1,ki - 1
                       work(k + iv*n) = czero
                    end do
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == nb) .or. (ki == n)) then
                       call la_zgemm('N','N',n,iv,n - ki + iv,cone,vl(1,ki - iv + 1),ldvl, &
                                 work(ki - iv + 1 + (1)*n),n,czero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          ii = la_izamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_zdscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_zlacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ztrevc3
#ifdef LA_WITH_XDP
     !> YTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by YHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_ytrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,rwork,lrwork,info)
        use la_constants_xdp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,lrwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki,iv,maxwrk,nb
           real(xdp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(xdp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           nb = la_ilaenv(1,'YTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           rwork(1) = n
           lquery = (lwork == -1 .or. lrwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -14
           else if (lrwork < max(1,n) .and. .not. lquery) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('YTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_ylaset('F',n,1 + 2*nb,czero,czero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_xlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_xlabad(unfl,ovfl)
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_xyasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=nb=1;
              ! blocked     version starts with iv=nb, goes down to 1.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = nb
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex right eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k + iv*n) = -t(k,ki)
                 end do
                 ! solve upper triangular system:
                 ! [ t(1:ki-1,1:ki-1) - t(ki,ki) ]*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_ylatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vr and normalize.
                    call la_ycopy(ki,work(1 + iv*n),1,vr(1,is),1)
                    ii = la_iyamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_yxscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki > 1) call la_ygemv('N',n,ki - 1,cone,vr,ldvr,work(1 + iv*n),1, &
                               cmplx(scale,KIND=xdp),vr(1,ki),1)
                    ii = la_iyamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_yxscal(n,remax,vr(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out below vector
                    do k = ki + 1,n
                       work(k + iv*n) = czero
                    end do
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == 1) .or. (ki == 1)) then
                       call la_ygemm('N','N',n,nb - iv + 1,ki + nb - iv,cone,vr,ldvr,work(1 + &
                                 (iv)*n),n,czero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          ii = la_iyamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_yxscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_ylacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = 1
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex left eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k + iv*n) = -conjg(t(ki,k))
                 end do
                 ! solve conjugate-transposed triangular system:
                 ! [ t(ki+1:n,ki+1:n) - t(ki,ki) ]**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_ylatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vl and normalize.
                    call la_ycopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                    ii = la_iyamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_yxscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki < n) call la_ygemv('N',n,n - ki,cone,vl(1,ki + 1),ldvl,work(ki + &
                              1 + iv*n),1,cmplx(scale,KIND=xdp),vl(1,ki),1)
                    ii = la_iyamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_yxscal(n,remax,vl(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out above vector
                    ! could go from ki-nv+1 to ki-1
                    do k = 1,ki - 1
                       work(k + iv*n) = czero
                    end do
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == nb) .or. (ki == n)) then
                       call la_ygemm('N','N',n,iv,n - ki + iv,cone,vl(1,ki - iv + 1),ldvl, &
                                 work(ki - iv + 1 + (1)*n),n,czero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          ii = la_iyamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_yxscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_ylacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_ytrevc3
#endif
#ifdef LA_WITH_QP
     !> WTREVC3: computes some or all of the right and/or left eigenvectors of
     !> a complex upper triangular matrix T.
     !> Matrices of this type are produced by the Schur factorization of
     !> a complex general matrix:  A = Q*T*Q**H, as computed by WHSEQR.
     !> The right eigenvector x and the left eigenvector y of T corresponding
     !> to an eigenvalue w are defined by:
     !> T*x = w*x,     (y**H)*T = w*(y**H)
     !> where y**H denotes the conjugate transpose of the vector y.
     !> The eigenvalues are not input to this routine, but are read directly
     !> from the diagonal of T.
     !> This routine returns the matrices X and/or Y of right and left
     !> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
     !> input matrix. If Q is the unitary factor that reduces a matrix A to
     !> Schur form T, then Q*X and Q*Y are the matrices of right and left
     !> eigenvectors of A.
     !> This uses a Level 3 BLAS version of the back transformation.

     pure subroutine la_wtrevc3(side,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,mm,m, &
               work,lwork,rwork,lrwork,info)
        use la_constants_qp,only:zero,one,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,lwork,lrwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmin = 8
           integer(ilp),parameter :: nbmax = 128

           ! Local Scalars
           logical(lk) :: allv,bothv,leftv,lquery,over,rightv,somev
           integer(ilp) :: i,ii,is,j,k,ki,iv,maxwrk,nb
           real(qp) :: ovfl,remax,scale,smin,smlnum,ulp,unfl
           complex(qp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           allv = la_lsame(howmny,'A')
           over = la_lsame(howmny,'B')
           somev = la_lsame(howmny,'S')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           if (somev) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           nb = la_ilaenv(1,'WTREVC',side//howmny,n,-1,-1,-1)
           maxwrk = n + 2*n*nb
           work(1) = maxwrk
           rwork(1) = n
           lquery = (lwork == -1 .or. lrwork == -1)
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. allv .and. .not. over .and. .not. somev) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -11
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -14
           else if (lrwork < max(1,n) .and. .not. lquery) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('WTREVC3',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! use blocked version of back-transformation if sufficient workspace.
           ! zero-out the workspace to avoid potential nan propagation.
           if (over .and. lwork >= n + 2*n*nbmin) then
              nb = (lwork - n)/(2*n)
              nb = min(nb,nbmax)
              call la_wlaset('F',n,1 + 2*nb,czero,czero,work,n)
           else
              nb = 1
           end if
           ! set the constants to control overflow.
           unfl = la_qlamch('SAFE MINIMUM')
           ovfl = one/unfl
           call la_qlabad(unfl,ovfl)
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ! store the diagonal elements of t in working array work.
           do i = 1,n
              work(i) = t(i,i)
           end do
           ! compute 1-norm of each column of strictly upper triangular
           ! part of t to control overflow in triangular solver.
           rwork(1) = zero
           do j = 2,n
              rwork(j) = la_qwasum(j - 1,t(1,j),1)
           end do
           if (rightv) then
              ! ============================================================
              ! compute right eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=nb=1;
              ! blocked     version starts with iv=nb, goes down to 1.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = nb
              is = m
              loop_80: do ki = n,1,-1
                 if (somev) then
                    if (.not. select(ki)) cycle loop_80
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex right eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = 1,ki - 1
                    work(k + iv*n) = -t(k,ki)
                 end do
                 ! solve upper triangular system:
                 ! [ t(1:ki-1,1:ki-1) - t(ki,ki) ]*x = scale*work.
                 do k = 1,ki - 1
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki > 1) then
                    call la_wlatrs('UPPER','NO TRANSPOSE','NON-UNIT','Y',ki - 1,t,ldt, &
                              work(1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vr and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vr and normalize.
                    call la_wcopy(ki,work(1 + iv*n),1,vr(1,is),1)
                    ii = la_iwamax(ki,vr(1,is),1)
                    remax = one/cabs1(vr(ii,is))
                    call la_wqscal(ki,remax,vr(1,is),1)
                    do k = ki + 1,n
                       vr(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki > 1) call la_wgemv('N',n,ki - 1,cone,vr,ldvr,work(1 + iv*n),1, &
                               cmplx(scale,KIND=qp),vr(1,ki),1)
                    ii = la_iwamax(n,vr(1,ki),1)
                    remax = one/cabs1(vr(ii,ki))
                    call la_wqscal(n,remax,vr(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out below vector
                    do k = ki + 1,n
                       work(k + iv*n) = czero
                    end do
                    ! columns iv:nb of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == 1) .or. (ki == 1)) then
                       call la_wgemm('N','N',n,nb - iv + 1,ki + nb - iv,cone,vr,ldvr,work(1 + &
                                 (iv)*n),n,czero,work(1 + (nb + iv)*n),n)
                       ! normalize vectors
                       do k = iv,nb
                          ii = la_iwamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_wqscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_wlacpy('F',n,nb - iv + 1,work(1 + (nb + iv)*n),n,vr(1,ki), &
                                 ldvr)
                       iv = nb
                    else
                       iv = iv - 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = 1,ki - 1
                    t(k,k) = work(k)
                 end do
                 is = is - 1
              end do loop_80
           end if
           if (leftv) then
              ! ============================================================
              ! compute left eigenvectors.
              ! iv is index of column in current block.
              ! non-blocked version always uses iv=1;
              ! blocked     version starts with iv=1, goes up to nb.
              ! (note the "0-th" column is used to store the original diagonal.)
              iv = 1
              is = 1
              loop_130: do ki = 1,n
                 if (somev) then
                    if (.not. select(ki)) cycle loop_130
                 end if
                 smin = max(ulp*(cabs1(t(ki,ki))),smlnum)
                 ! --------------------------------------------------------
                 ! complex left eigenvector
                 work(ki + iv*n) = cone
                 ! form right-hand side.
                 do k = ki + 1,n
                    work(k + iv*n) = -conjg(t(ki,k))
                 end do
                 ! solve conjugate-transposed triangular system:
                 ! [ t(ki+1:n,ki+1:n) - t(ki,ki) ]**h * x = scale*work.
                 do k = ki + 1,n
                    t(k,k) = t(k,k) - t(ki,ki)
                    if (cabs1(t(k,k)) < smin) t(k,k) = smin
                 end do
                 if (ki < n) then
                    call la_wlatrs('UPPER','CONJUGATE TRANSPOSE','NON-UNIT','Y',n - ki,t( &
                              ki + 1,ki + 1),ldt,work(ki + 1 + iv*n),scale,rwork,info)
                    work(ki + iv*n) = scale
                 end if
                 ! copy the vector x or q*x to vl and normalize.
                 if (.not. over) then
                    ! ------------------------------
                    ! no back-transform: copy x to vl and normalize.
                    call la_wcopy(n - ki + 1,work(ki + iv*n),1,vl(ki,is),1)
                    ii = la_iwamax(n - ki + 1,vl(ki,is),1) + ki - 1
                    remax = one/cabs1(vl(ii,is))
                    call la_wqscal(n - ki + 1,remax,vl(ki,is),1)
                    do k = 1,ki - 1
                       vl(k,is) = czero
                    end do
                 else if (nb == 1) then
                    ! ------------------------------
                    ! version 1: back-transform each vector with gemv, q*x.
                    if (ki < n) call la_wgemv('N',n,n - ki,cone,vl(1,ki + 1),ldvl,work(ki + &
                              1 + iv*n),1,cmplx(scale,KIND=qp),vl(1,ki),1)
                    ii = la_iwamax(n,vl(1,ki),1)
                    remax = one/cabs1(vl(ii,ki))
                    call la_wqscal(n,remax,vl(1,ki),1)
                 else
                    ! ------------------------------
                    ! version 2: back-transform block of vectors with gemm
                    ! zero out above vector
                    ! could go from ki-nv+1 to ki-1
                    do k = 1,ki - 1
                       work(k + iv*n) = czero
                    end do
                    ! columns 1:iv of work are valid vectors.
                    ! when the number of vectors stored reaches nb,
                    ! or if this was last vector, do the gemm
                    if ((iv == nb) .or. (ki == n)) then
                       call la_wgemm('N','N',n,iv,n - ki + iv,cone,vl(1,ki - iv + 1),ldvl, &
                                 work(ki - iv + 1 + (1)*n),n,czero,work(1 + (nb + 1)*n),n)
                       ! normalize vectors
                       do k = 1,iv
                          ii = la_iwamax(n,work(1 + (nb + k)*n),1)
                          remax = one/cabs1(work(ii + (nb + k)*n))
                          call la_wqscal(n,remax,work(1 + (nb + k)*n),1)
                       end do
                       call la_wlacpy('F',n,iv,work(1 + (nb + 1)*n),n,vl(1,ki - iv + 1), &
                                 ldvl)
                       iv = 1
                    else
                       iv = iv + 1
                    end if
                 end if
                 ! restore the original diagonal elements of t.
                 do k = ki + 1,n
                    t(k,k) = work(k)
                 end do
                 is = is + 1
              end do loop_130
           end if
           return
     end subroutine la_wtrevc3
#endif

     !> CTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a complex upper triangular
     !> matrix T (or of any matrix Q*T*Q**H with Q unitary).

     pure subroutine la_ctrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm, &
                m,work,ldwork,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(sp),intent(out) :: rwork(*),s(*),sep(*)
           complex(sp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(sp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: somcon,wantbh,wants,wantsp
           character :: normin
           integer(ilp) :: i,ierr,ix,j,k,kase,ks
           real(sp) :: bignum,eps,est,lnrm,rnrm,scale,smlnum,xnorm
           complex(sp) :: cdum,prod
           ! Local Arrays
           integer(ilp) :: isave(3)
           complex(sp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           ! set m to the number of eigenpairs for which condition numbers are
           ! to be computed.
           if (somcon) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -13
           else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('CTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ks = 1
           loop_50: do k = 1,n
              if (somcon) then
                 if (.not. select(k)) cycle loop_50
              end if
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 prod = la_cdotc(n,vr(1,ks),1,vl(1,ks),1)
                 rnrm = la_scnrm2(n,vr(1,ks),1)
                 lnrm = la_scnrm2(n,vl(1,ks),1)
                 s(ks) = abs(prod)/(rnrm*lnrm)
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the k-th
                 ! diagonal element to the (1,1) position.
                 call la_clacpy('FULL',n,n,t,ldt,work,ldwork)
                 call la_ctrexc('NO Q',n,work,ldwork,dummy,1,k,1,ierr)
                 ! form  c = t22 - lambda*i in work(2:n,2:n).
                 do i = 2,n
                    work(i,i) = work(i,i) - work(1,1)
                 end do
                 ! estimate a lower bound for the 1-norm of inv(c**h). the 1st
                 ! and (n+1)th columns of work are used to store work vectors.
                 sep(ks) = zero
                 est = zero
                 kase = 0
                 normin = 'N'
                 30 continue
                 call la_clacn2(n - 1,work(1,n + 1),work,est,kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve c**h*x = scale*b
                       call la_clatrs('UPPER','CONJUGATE TRANSPOSE','NONUNIT',normin,n - 1, &
                                 work(2,2),ldwork,work,scale,rwork,ierr)
                    else
                       ! solve c*x = scale*b
                       call la_clatrs('UPPER','NO TRANSPOSE','NONUNIT',normin,n - 1,work( &
                                 2,2),ldwork,work,scale,rwork,ierr)
                    end if
                    normin = 'Y'
                    if (scale /= one) then
                       ! multiply by 1/scale if doing so will not cause
                       ! overflow.
                       ix = la_icamax(n - 1,work,1)
                       xnorm = cabs1(work(ix,1))
                       if (scale < xnorm*smlnum .or. scale == zero) go to 40
                       call la_csrscl(n,scale,work,1)
                    end if
                    go to 30
                 end if
                 sep(ks) = one/max(est,smlnum)
              end if
              40 continue
              ks = ks + 1
           end do loop_50
           return
     end subroutine la_ctrsna
     !> ZTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a complex upper triangular
     !> matrix T (or of any matrix Q*T*Q**H with Q unitary).

     pure subroutine la_ztrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm, &
                m,work,ldwork,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(dp),intent(out) :: rwork(*),s(*),sep(*)
           complex(dp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(dp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: somcon,wantbh,wants,wantsp
           character :: normin
           integer(ilp) :: i,ierr,ix,j,k,kase,ks
           real(dp) :: bignum,eps,est,lnrm,rnrm,scale,smlnum,xnorm
           complex(dp) :: cdum,prod
           ! Local Arrays
           integer(ilp) :: isave(3)
           complex(dp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           ! set m to the number of eigenpairs for which condition numbers are
           ! to be computed.
           if (somcon) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -13
           else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('ZTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ks = 1
           loop_50: do k = 1,n
              if (somcon) then
                 if (.not. select(k)) cycle loop_50
              end if
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 prod = la_zdotc(n,vr(1,ks),1,vl(1,ks),1)
                 rnrm = la_dznrm2(n,vr(1,ks),1)
                 lnrm = la_dznrm2(n,vl(1,ks),1)
                 s(ks) = abs(prod)/(rnrm*lnrm)
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the k-th
                 ! diagonal element to the (1,1) position.
                 call la_zlacpy('FULL',n,n,t,ldt,work,ldwork)
                 call la_ztrexc('NO Q',n,work,ldwork,dummy,1,k,1,ierr)
                 ! form  c = t22 - lambda*i in work(2:n,2:n).
                 do i = 2,n
                    work(i,i) = work(i,i) - work(1,1)
                 end do
                 ! estimate a lower bound for the 1-norm of inv(c**h). the 1st
                 ! and (n+1)th columns of work are used to store work vectors.
                 sep(ks) = zero
                 est = zero
                 kase = 0
                 normin = 'N'
                 30 continue
                 call la_zlacn2(n - 1,work(1,n + 1),work,est,kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve c**h*x = scale*b
                       call la_zlatrs('UPPER','CONJUGATE TRANSPOSE','NONUNIT',normin,n - 1, &
                                 work(2,2),ldwork,work,scale,rwork,ierr)
                    else
                       ! solve c*x = scale*b
                       call la_zlatrs('UPPER','NO TRANSPOSE','NONUNIT',normin,n - 1,work( &
                                 2,2),ldwork,work,scale,rwork,ierr)
                    end if
                    normin = 'Y'
                    if (scale /= one) then
                       ! multiply by 1/scale if doing so will not cause
                       ! overflow.
                       ix = la_izamax(n - 1,work,1)
                       xnorm = cabs1(work(ix,1))
                       if (scale < xnorm*smlnum .or. scale == zero) go to 40
                       call la_zdrscl(n,scale,work,1)
                    end if
                    go to 30
                 end if
                 sep(ks) = one/max(est,smlnum)
              end if
              40 continue
              ks = ks + 1
           end do loop_50
           return
     end subroutine la_ztrsna
#ifdef LA_WITH_XDP
     !> YTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a complex upper triangular
     !> matrix T (or of any matrix Q*T*Q**H with Q unitary).

     pure subroutine la_ytrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm, &
                m,work,ldwork,rwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(xdp),intent(out) :: rwork(*),s(*),sep(*)
           complex(xdp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(xdp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: somcon,wantbh,wants,wantsp
           character :: normin
           integer(ilp) :: i,ierr,ix,j,k,kase,ks
           real(xdp) :: bignum,eps,est,lnrm,rnrm,scale,smlnum,xnorm
           complex(xdp) :: cdum,prod
           ! Local Arrays
           integer(ilp) :: isave(3)
           complex(xdp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           ! set m to the number of eigenpairs for which condition numbers are
           ! to be computed.
           if (somcon) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -13
           else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('YTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ks = 1
           loop_50: do k = 1,n
              if (somcon) then
                 if (.not. select(k)) cycle loop_50
              end if
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 prod = la_ydotc(n,vr(1,ks),1,vl(1,ks),1)
                 rnrm = la_xynrm2(n,vr(1,ks),1)
                 lnrm = la_xynrm2(n,vl(1,ks),1)
                 s(ks) = abs(prod)/(rnrm*lnrm)
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the k-th
                 ! diagonal element to the (1,1) position.
                 call la_ylacpy('FULL',n,n,t,ldt,work,ldwork)
                 call la_ytrexc('NO Q',n,work,ldwork,dummy,1,k,1,ierr)
                 ! form  c = t22 - lambda*i in work(2:n,2:n).
                 do i = 2,n
                    work(i,i) = work(i,i) - work(1,1)
                 end do
                 ! estimate a lower bound for the 1-norm of inv(c**h). the 1st
                 ! and (n+1)th columns of work are used to store work vectors.
                 sep(ks) = zero
                 est = zero
                 kase = 0
                 normin = 'N'
                 30 continue
                 call la_ylacn2(n - 1,work(1,n + 1),work,est,kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve c**h*x = scale*b
                       call la_ylatrs('UPPER','CONJUGATE TRANSPOSE','NONUNIT',normin,n - 1, &
                                 work(2,2),ldwork,work,scale,rwork,ierr)
                    else
                       ! solve c*x = scale*b
                       call la_ylatrs('UPPER','NO TRANSPOSE','NONUNIT',normin,n - 1,work( &
                                 2,2),ldwork,work,scale,rwork,ierr)
                    end if
                    normin = 'Y'
                    if (scale /= one) then
                       ! multiply by 1/scale if doing so will not cause
                       ! overflow.
                       ix = la_iyamax(n - 1,work,1)
                       xnorm = cabs1(work(ix,1))
                       if (scale < xnorm*smlnum .or. scale == zero) go to 40
                       call la_yxrscl(n,scale,work,1)
                    end if
                    go to 30
                 end if
                 sep(ks) = one/max(est,smlnum)
              end if
              40 continue
              ks = ks + 1
           end do loop_50
           return
     end subroutine la_ytrsna
#endif
#ifdef LA_WITH_QP
     !> WTRSNA: estimates reciprocal condition numbers for specified
     !> eigenvalues and/or right eigenvectors of a complex upper triangular
     !> matrix T (or of any matrix Q*T*Q**H with Q unitary).

     pure subroutine la_wtrsna(job,howmny,select,n,t,ldt,vl,ldvl,vr,ldvr,s,sep,mm, &
                m,work,ldwork,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: howmny,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldt,ldvl,ldvr,ldwork,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           real(qp),intent(out) :: rwork(*),s(*),sep(*)
           complex(qp),intent(in) :: t(ldt,*),vl(ldvl,*),vr(ldvr,*)
           complex(qp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: somcon,wantbh,wants,wantsp
           character :: normin
           integer(ilp) :: i,ierr,ix,j,k,kase,ks
           real(qp) :: bignum,eps,est,lnrm,rnrm,scale,smlnum,xnorm
           complex(qp) :: cdum,prod
           ! Local Arrays
           integer(ilp) :: isave(3)
           complex(qp) :: dummy(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           somcon = la_lsame(howmny,'S')
           ! set m to the number of eigenpairs for which condition numbers are
           ! to be computed.
           if (somcon) then
              m = 0
              do j = 1,n
                 if (select(j)) m = m + 1
              end do
           else
              m = n
           end if
           info = 0
           if (.not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(howmny,'A') .and. .not. somcon) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldvl < 1 .or. (wants .and. ldvl < n)) then
              info = -8
           else if (ldvr < 1 .or. (wants .and. ldvr < n)) then
              info = -10
           else if (mm < m) then
              info = -13
           else if (ldwork < 1 .or. (wantsp .and. ldwork < n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('WTRSNA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (somcon) then
                 if (.not. select(1)) return
              end if
              if (wants) s(1) = one
              if (wantsp) sep(1) = abs(t(1,1))
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ks = 1
           loop_50: do k = 1,n
              if (somcon) then
                 if (.not. select(k)) cycle loop_50
              end if
              if (wants) then
                 ! compute the reciprocal condition number of the k-th
                 ! eigenvalue.
                 prod = la_wdotc(n,vr(1,ks),1,vl(1,ks),1)
                 rnrm = la_qwnrm2(n,vr(1,ks),1)
                 lnrm = la_qwnrm2(n,vl(1,ks),1)
                 s(ks) = abs(prod)/(rnrm*lnrm)
              end if
              if (wantsp) then
                 ! estimate the reciprocal condition number of the k-th
                 ! eigenvector.
                 ! copy the matrix t to the array work and swap the k-th
                 ! diagonal element to the (1,1) position.
                 call la_wlacpy('FULL',n,n,t,ldt,work,ldwork)
                 call la_wtrexc('NO Q',n,work,ldwork,dummy,1,k,1,ierr)
                 ! form  c = t22 - lambda*i in work(2:n,2:n).
                 do i = 2,n
                    work(i,i) = work(i,i) - work(1,1)
                 end do
                 ! estimate a lower bound for the 1-norm of inv(c**h). the 1st
                 ! and (n+1)th columns of work are used to store work vectors.
                 sep(ks) = zero
                 est = zero
                 kase = 0
                 normin = 'N'
                 30 continue
                 call la_wlacn2(n - 1,work(1,n + 1),work,est,kase,isave)
                 if (kase /= 0) then
                    if (kase == 1) then
                       ! solve c**h*x = scale*b
                       call la_wlatrs('UPPER','CONJUGATE TRANSPOSE','NONUNIT',normin,n - 1, &
                                 work(2,2),ldwork,work,scale,rwork,ierr)
                    else
                       ! solve c*x = scale*b
                       call la_wlatrs('UPPER','NO TRANSPOSE','NONUNIT',normin,n - 1,work( &
                                 2,2),ldwork,work,scale,rwork,ierr)
                    end if
                    normin = 'Y'
                    if (scale /= one) then
                       ! multiply by 1/scale if doing so will not cause
                       ! overflow.
                       ix = la_iwamax(n - 1,work,1)
                       xnorm = cabs1(work(ix,1))
                       if (scale < xnorm*smlnum .or. scale == zero) go to 40
                       call la_wqrscl(n,scale,work,1)
                    end if
                    go to 30
                 end if
                 sep(ks) = one/max(est,smlnum)
              end if
              40 continue
              ks = ks + 1
           end do loop_50
           return
     end subroutine la_wtrsna
#endif

     !> CLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue W of a complex upper Hessenberg
     !> matrix H.

     pure subroutine la_claein(rightv,noinit,n,h,ldh,w,v,b,ldb,rwork,eps3,smlnum, &
               info)
        use la_constants_sp,only:one,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(sp),intent(in) :: eps3,smlnum
           complex(sp),intent(in) :: w
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(out) :: b(ldb,*)
           complex(sp),intent(in) :: h(ldh,*)
           complex(sp),intent(inout) :: v(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: tenth = 1.0e-1_sp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,ierr,its,j
           real(sp) :: growto,nrmsml,rootn,rtemp,scale,vnorm
           complex(sp) :: cdum,ei,ej,temp,x
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=sp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - w*i (except that the subdiagonal elements are not
           ! stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - w
           end do
           if (noinit) then
              ! initialize v.
              do i = 1,n
                 v(i) = eps3
              end do
           else
              ! scale supplied initial vector.
              vnorm = la_scnrm2(n,v,1)
              call la_csscal(n, (eps3*rootn)/max(vnorm,nrmsml),v,1)
           end if
           if (rightv) then
              ! lu decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do i = 1,n - 1
                 ei = h(i + 1,i)
                 if (cabs1(b(i,i)) < cabs1(ei)) then
                    ! interchange rows and eliminate.
                    x = la_cladiv(b(i,i),ei)
                    b(i,i) = ei
                    do j = i + 1,n
                       temp = b(i + 1,j)
                       b(i + 1,j) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(i,i) == czero) b(i,i) = eps3
                    x = la_cladiv(ei,b(i,i))
                    if (x /= czero) then
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(n,n) == czero) b(n,n) = eps3
              trans = 'N'
           else
              ! ul decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do j = n,2,-1
                 ej = h(j,j - 1)
                 if (cabs1(b(j,j)) < cabs1(ej)) then
                    ! interchange columns and eliminate.
                    x = la_cladiv(b(j,j),ej)
                    b(j,j) = ej
                    do i = 1,j - 1
                       temp = b(i,j - 1)
                       b(i,j - 1) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(j,j) == czero) b(j,j) = eps3
                    x = la_cladiv(ej,b(j,j))
                    if (x /= czero) then
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(1,1) == czero) b(1,1) = eps3
              trans = 'C'
           end if
           normin = 'N'
           do its = 1,n
              ! solve u*x = scale*v for a right eigenvector
                ! or u**h *x = scale*v for a left eigenvector,
              ! overwriting x on v.
              call la_clatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,v,scale,rwork, &
                        ierr)
              normin = 'Y'
              ! test for sufficient growth in the norm of v.
              vnorm = la_scasum(n,v,1)
              if (vnorm >= growto*scale) go to 120
              ! choose new orthogonal starting vector and try again.
              rtemp = eps3/(rootn + one)
              v(1) = eps3
              do i = 2,n
                 v(i) = rtemp
              end do
              v(n - its + 1) = v(n - its + 1) - eps3*rootn
           end do
           ! failure to find eigenvector in n iterations.
           info = 1
           120 continue
           ! normalize eigenvector.
           i = la_icamax(n,v,1)
           call la_csscal(n,one/cabs1(v(i)),v,1)
           return
     end subroutine la_claein
     !> ZLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue W of a complex upper Hessenberg
     !> matrix H.

     pure subroutine la_zlaein(rightv,noinit,n,h,ldh,w,v,b,ldb,rwork,eps3,smlnum, &
               info)
        use la_constants_dp,only:one,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(dp),intent(in) :: eps3,smlnum
           complex(dp),intent(in) :: w
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: b(ldb,*)
           complex(dp),intent(in) :: h(ldh,*)
           complex(dp),intent(inout) :: v(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: tenth = 1.0e-1_dp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,ierr,its,j
           real(dp) :: growto,nrmsml,rootn,rtemp,scale,vnorm
           complex(dp) :: cdum,ei,ej,temp,x
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=dp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - w*i (except that the subdiagonal elements are not
           ! stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - w
           end do
           if (noinit) then
              ! initialize v.
              do i = 1,n
                 v(i) = eps3
              end do
           else
              ! scale supplied initial vector.
              vnorm = la_dznrm2(n,v,1)
              call la_zdscal(n, (eps3*rootn)/max(vnorm,nrmsml),v,1)
           end if
           if (rightv) then
              ! lu decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do i = 1,n - 1
                 ei = h(i + 1,i)
                 if (cabs1(b(i,i)) < cabs1(ei)) then
                    ! interchange rows and eliminate.
                    x = la_zladiv(b(i,i),ei)
                    b(i,i) = ei
                    do j = i + 1,n
                       temp = b(i + 1,j)
                       b(i + 1,j) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(i,i) == czero) b(i,i) = eps3
                    x = la_zladiv(ei,b(i,i))
                    if (x /= czero) then
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(n,n) == czero) b(n,n) = eps3
              trans = 'N'
           else
              ! ul decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do j = n,2,-1
                 ej = h(j,j - 1)
                 if (cabs1(b(j,j)) < cabs1(ej)) then
                    ! interchange columns and eliminate.
                    x = la_zladiv(b(j,j),ej)
                    b(j,j) = ej
                    do i = 1,j - 1
                       temp = b(i,j - 1)
                       b(i,j - 1) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(j,j) == czero) b(j,j) = eps3
                    x = la_zladiv(ej,b(j,j))
                    if (x /= czero) then
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(1,1) == czero) b(1,1) = eps3
              trans = 'C'
           end if
           normin = 'N'
           do its = 1,n
              ! solve u*x = scale*v for a right eigenvector
                ! or u**h *x = scale*v for a left eigenvector,
              ! overwriting x on v.
              call la_zlatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,v,scale,rwork, &
                        ierr)
              normin = 'Y'
              ! test for sufficient growth in the norm of v.
              vnorm = la_dzasum(n,v,1)
              if (vnorm >= growto*scale) go to 120
              ! choose new orthogonal starting vector and try again.
              rtemp = eps3/(rootn + one)
              v(1) = eps3
              do i = 2,n
                 v(i) = rtemp
              end do
              v(n - its + 1) = v(n - its + 1) - eps3*rootn
           end do
           ! failure to find eigenvector in n iterations.
           info = 1
           120 continue
           ! normalize eigenvector.
           i = la_izamax(n,v,1)
           call la_zdscal(n,one/cabs1(v(i)),v,1)
           return
     end subroutine la_zlaein
#ifdef LA_WITH_XDP
     !> YLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue W of a complex upper Hessenberg
     !> matrix H.

     pure subroutine la_ylaein(rightv,noinit,n,h,ldh,w,v,b,ldb,rwork,eps3,smlnum, &
               info)
        use la_constants_xdp,only:one,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(xdp),intent(in) :: eps3,smlnum
           complex(xdp),intent(in) :: w
           ! Array Arguments
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(out) :: b(ldb,*)
           complex(xdp),intent(in) :: h(ldh,*)
           complex(xdp),intent(inout) :: v(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: tenth = 1.0e-1_xdp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,ierr,its,j
           real(xdp) :: growto,nrmsml,rootn,rtemp,scale,vnorm
           complex(xdp) :: cdum,ei,ej,temp,x
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=xdp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - w*i (except that the subdiagonal elements are not
           ! stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - w
           end do
           if (noinit) then
              ! initialize v.
              do i = 1,n
                 v(i) = eps3
              end do
           else
              ! scale supplied initial vector.
              vnorm = la_xynrm2(n,v,1)
              call la_yxscal(n, (eps3*rootn)/max(vnorm,nrmsml),v,1)
           end if
           if (rightv) then
              ! lu decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do i = 1,n - 1
                 ei = h(i + 1,i)
                 if (cabs1(b(i,i)) < cabs1(ei)) then
                    ! interchange rows and eliminate.
                    x = la_yladiv(b(i,i),ei)
                    b(i,i) = ei
                    do j = i + 1,n
                       temp = b(i + 1,j)
                       b(i + 1,j) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(i,i) == czero) b(i,i) = eps3
                    x = la_yladiv(ei,b(i,i))
                    if (x /= czero) then
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(n,n) == czero) b(n,n) = eps3
              trans = 'N'
           else
              ! ul decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do j = n,2,-1
                 ej = h(j,j - 1)
                 if (cabs1(b(j,j)) < cabs1(ej)) then
                    ! interchange columns and eliminate.
                    x = la_yladiv(b(j,j),ej)
                    b(j,j) = ej
                    do i = 1,j - 1
                       temp = b(i,j - 1)
                       b(i,j - 1) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(j,j) == czero) b(j,j) = eps3
                    x = la_yladiv(ej,b(j,j))
                    if (x /= czero) then
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(1,1) == czero) b(1,1) = eps3
              trans = 'C'
           end if
           normin = 'N'
           do its = 1,n
              ! solve u*x = scale*v for a right eigenvector
                ! or u**h *x = scale*v for a left eigenvector,
              ! overwriting x on v.
              call la_ylatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,v,scale,rwork, &
                        ierr)
              normin = 'Y'
              ! test for sufficient growth in the norm of v.
              vnorm = la_xyasum(n,v,1)
              if (vnorm >= growto*scale) go to 120
              ! choose new orthogonal starting vector and try again.
              rtemp = eps3/(rootn + one)
              v(1) = eps3
              do i = 2,n
                 v(i) = rtemp
              end do
              v(n - its + 1) = v(n - its + 1) - eps3*rootn
           end do
           ! failure to find eigenvector in n iterations.
           info = 1
           120 continue
           ! normalize eigenvector.
           i = la_iyamax(n,v,1)
           call la_yxscal(n,one/cabs1(v(i)),v,1)
           return
     end subroutine la_ylaein
#endif
#ifdef LA_WITH_QP
     !> WLAEIN: uses inverse iteration to find a right or left eigenvector
     !> corresponding to the eigenvalue W of a complex upper Hessenberg
     !> matrix H.

     pure subroutine la_wlaein(rightv,noinit,n,h,ldh,w,v,b,ldb,rwork,eps3,smlnum, &
               info)
        use la_constants_qp,only:one,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: noinit,rightv
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldh,n
           real(qp),intent(in) :: eps3,smlnum
           complex(qp),intent(in) :: w
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(out) :: b(ldb,*)
           complex(qp),intent(in) :: h(ldh,*)
           complex(qp),intent(inout) :: v(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: tenth = 1.0e-1_qp

           ! Local Scalars
           character :: normin,trans
           integer(ilp) :: i,ierr,its,j
           real(qp) :: growto,nrmsml,rootn,rtemp,scale,vnorm
           complex(qp) :: cdum,ei,ej,temp,x
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           info = 0
           ! growto is the threshold used in the acceptance test for an
           ! eigenvector.
           rootn = sqrt(real(n,KIND=qp))
           growto = tenth/rootn
           nrmsml = max(one,eps3*rootn)*smlnum
           ! form b = h - w*i (except that the subdiagonal elements are not
           ! stored).
           do j = 1,n
              do i = 1,j - 1
                 b(i,j) = h(i,j)
              end do
              b(j,j) = h(j,j) - w
           end do
           if (noinit) then
              ! initialize v.
              do i = 1,n
                 v(i) = eps3
              end do
           else
              ! scale supplied initial vector.
              vnorm = la_qwnrm2(n,v,1)
              call la_wqscal(n, (eps3*rootn)/max(vnorm,nrmsml),v,1)
           end if
           if (rightv) then
              ! lu decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do i = 1,n - 1
                 ei = h(i + 1,i)
                 if (cabs1(b(i,i)) < cabs1(ei)) then
                    ! interchange rows and eliminate.
                    x = la_wladiv(b(i,i),ei)
                    b(i,i) = ei
                    do j = i + 1,n
                       temp = b(i + 1,j)
                       b(i + 1,j) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(i,i) == czero) b(i,i) = eps3
                    x = la_wladiv(ei,b(i,i))
                    if (x /= czero) then
                       do j = i + 1,n
                          b(i + 1,j) = b(i + 1,j) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(n,n) == czero) b(n,n) = eps3
              trans = 'N'
           else
              ! ul decomposition with partial pivoting of b, replacing czero
              ! pivots by eps3.
              do j = n,2,-1
                 ej = h(j,j - 1)
                 if (cabs1(b(j,j)) < cabs1(ej)) then
                    ! interchange columns and eliminate.
                    x = la_wladiv(b(j,j),ej)
                    b(j,j) = ej
                    do i = 1,j - 1
                       temp = b(i,j - 1)
                       b(i,j - 1) = b(i,j) - x*temp
                       b(i,j) = temp
                    end do
                 else
                    ! eliminate without interchange.
                    if (b(j,j) == czero) b(j,j) = eps3
                    x = la_wladiv(ej,b(j,j))
                    if (x /= czero) then
                       do i = 1,j - 1
                          b(i,j - 1) = b(i,j - 1) - x*b(i,j)
                       end do
                    end if
                 end if
              end do
              if (b(1,1) == czero) b(1,1) = eps3
              trans = 'C'
           end if
           normin = 'N'
           do its = 1,n
              ! solve u*x = scale*v for a right eigenvector
                ! or u**h *x = scale*v for a left eigenvector,
              ! overwriting x on v.
              call la_wlatrs('UPPER',trans,'NONUNIT',normin,n,b,ldb,v,scale,rwork, &
                        ierr)
              normin = 'Y'
              ! test for sufficient growth in the norm of v.
              vnorm = la_qwasum(n,v,1)
              if (vnorm >= growto*scale) go to 120
              ! choose new orthogonal starting vector and try again.
              rtemp = eps3/(rootn + one)
              v(1) = eps3
              do i = 2,n
                 v(i) = rtemp
              end do
              v(n - its + 1) = v(n - its + 1) - eps3*rootn
           end do
           ! failure to find eigenvector in n iterations.
           info = 1
           120 continue
           ! normalize eigenvector.
           i = la_iwamax(n,v,1)
           call la_wqscal(n,one/cabs1(v(i)),v,1)
           return
     end subroutine la_wlaein
#endif

     !> CTRSYL: solves the complex Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**H, and A and B are both upper triangular. A is
     !> M-by-M and B is N-by-N; the right hand side C and the solution X are
     !> M-by-N; and scale is an output scale factor, set <= 1 to avoid
     !> overflow in X.

     subroutine la_ctrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_sp,only:one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: j,k,l
           real(sp) :: bignum,da11,db,eps,scaloc,sgn,smin,smlnum
           complex(sp) :: a11,suml,sumr,vec,x11
           ! Local Arrays
           real(sp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,min,real
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_slamch('P')
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=sp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_clange('M',m,m,a,lda,dum),eps*la_clange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                        l-1
                ! r(k,l) = sum [a(k,i)*x(i,l)] +isgn*sum [x(k,j)*b(j,l)].
                        ! i=k+1                      j=1
              loop_30: do l = 1,n
                 do k = m,1,-1
                    suml = la_cdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_cdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = a(k,k) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=sp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=sp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_cladiv(vec*cmplx(scaloc,KIND=sp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_csscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_30
           else if (.not. notrna .and. notrnb) then
              ! solve    a**h *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1                           l-1
                ! r(k,l) = sum [a**h(i,k)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                           j=1
              loop_60: do l = 1,n
                 do k = 1,m
                    suml = la_cdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_cdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = conjg(a(k,k)) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=sp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=sp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_cladiv(vec*cmplx(scaloc,KIND=sp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_csscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_60
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**h*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! upper-right corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! k-1
                 ! r(k,l) = sum [a**h(i,k)*x(i,l)] +
                          ! i=1
                                 ! n
                           ! isgn*sum [x(k,j)*b**h(l,j)].
                                ! j=l+1
              loop_90: do l = n,1,-1
                 do k = 1,m
                    suml = la_cdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_cdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = conjg(a(k,k) + sgn*b(l,l))
                    da11 = abs(real(a11,KIND=sp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=sp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_cladiv(vec*cmplx(scaloc,KIND=sp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_csscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_90
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                 ! a(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                          n
                ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b**h(l,j)]
                        ! i=k+1                      j=l+1
              loop_120: do l = n,1,-1
                 do k = m,1,-1
                    suml = la_cdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_cdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = a(k,k) + sgn*conjg(b(l,l))
                    da11 = abs(real(a11,KIND=sp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=sp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_cladiv(vec*cmplx(scaloc,KIND=sp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_csscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_120
           end if
           return
     end subroutine la_ctrsyl
     !> ZTRSYL: solves the complex Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**H, and A and B are both upper triangular. A is
     !> M-by-M and B is N-by-N; the right hand side C and the solution X are
     !> M-by-N; and scale is an output scale factor, set <= 1 to avoid
     !> overflow in X.

     subroutine la_ztrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_dp,only:one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: j,k,l
           real(dp) :: bignum,da11,db,eps,scaloc,sgn,smin,smlnum
           complex(dp) :: a11,suml,sumr,vec,x11
           ! Local Arrays
           real(dp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=dp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_zlange('M',m,m,a,lda,dum),eps*la_zlange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                        l-1
                ! r(k,l) = sum [a(k,i)*x(i,l)] +isgn*sum [x(k,j)*b(j,l)].
                        ! i=k+1                      j=1
              loop_30: do l = 1,n
                 do k = m,1,-1
                    suml = la_zdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_zdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = a(k,k) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=dp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=dp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_zladiv(vec*cmplx(scaloc,KIND=dp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_zdscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_30
           else if (.not. notrna .and. notrnb) then
              ! solve    a**h *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1                           l-1
                ! r(k,l) = sum [a**h(i,k)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                           j=1
              loop_60: do l = 1,n
                 do k = 1,m
                    suml = la_zdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_zdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = conjg(a(k,k)) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=dp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=dp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_zladiv(vec*cmplx(scaloc,KIND=dp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_zdscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_60
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**h*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! upper-right corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! k-1
                 ! r(k,l) = sum [a**h(i,k)*x(i,l)] +
                          ! i=1
                                 ! n
                           ! isgn*sum [x(k,j)*b**h(l,j)].
                                ! j=l+1
              loop_90: do l = n,1,-1
                 do k = 1,m
                    suml = la_zdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_zdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = conjg(a(k,k) + sgn*b(l,l))
                    da11 = abs(real(a11,KIND=dp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=dp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_zladiv(vec*cmplx(scaloc,KIND=dp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_zdscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_90
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                 ! a(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                          n
                ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b**h(l,j)]
                        ! i=k+1                      j=l+1
              loop_120: do l = n,1,-1
                 do k = m,1,-1
                    suml = la_zdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_zdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = a(k,k) + sgn*conjg(b(l,l))
                    da11 = abs(real(a11,KIND=dp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=dp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_zladiv(vec*cmplx(scaloc,KIND=dp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_zdscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_120
           end if
           return
     end subroutine la_ztrsyl
#ifdef LA_WITH_XDP
     !> YTRSYL: solves the complex Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**H, and A and B are both upper triangular. A is
     !> M-by-M and B is N-by-N; the right hand side C and the solution X are
     !> M-by-N; and scale is an output scale factor, set <= 1 to avoid
     !> overflow in X.

     subroutine la_ytrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_xdp,only:one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(xdp),intent(out) :: scale
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: j,k,l
           real(xdp) :: bignum,da11,db,eps,scaloc,sgn,smin,smlnum
           complex(xdp) :: a11,suml,sumr,vec,x11
           ! Local Arrays
           real(xdp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('YTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=xdp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_ylange('M',m,m,a,lda,dum),eps*la_ylange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                        l-1
                ! r(k,l) = sum [a(k,i)*x(i,l)] +isgn*sum [x(k,j)*b(j,l)].
                        ! i=k+1                      j=1
              loop_30: do l = 1,n
                 do k = m,1,-1
                    suml = la_ydotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_ydotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = a(k,k) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=xdp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=xdp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_yladiv(vec*cmplx(scaloc,KIND=xdp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_yxscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_30
           else if (.not. notrna .and. notrnb) then
              ! solve    a**h *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1                           l-1
                ! r(k,l) = sum [a**h(i,k)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                           j=1
              loop_60: do l = 1,n
                 do k = 1,m
                    suml = la_ydotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_ydotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = conjg(a(k,k)) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=xdp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=xdp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_yladiv(vec*cmplx(scaloc,KIND=xdp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_yxscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_60
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**h*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! upper-right corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! k-1
                 ! r(k,l) = sum [a**h(i,k)*x(i,l)] +
                          ! i=1
                                 ! n
                           ! isgn*sum [x(k,j)*b**h(l,j)].
                                ! j=l+1
              loop_90: do l = n,1,-1
                 do k = 1,m
                    suml = la_ydotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_ydotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = conjg(a(k,k) + sgn*b(l,l))
                    da11 = abs(real(a11,KIND=xdp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=xdp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_yladiv(vec*cmplx(scaloc,KIND=xdp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_yxscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_90
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                 ! a(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                          n
                ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b**h(l,j)]
                        ! i=k+1                      j=l+1
              loop_120: do l = n,1,-1
                 do k = m,1,-1
                    suml = la_ydotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_ydotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = a(k,k) + sgn*conjg(b(l,l))
                    da11 = abs(real(a11,KIND=xdp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=xdp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_yladiv(vec*cmplx(scaloc,KIND=xdp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_yxscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_120
           end if
           return
     end subroutine la_ytrsyl
#endif
#ifdef LA_WITH_QP
     !> WTRSYL: solves the complex Sylvester matrix equation:
     !> op(A)*X + X*op(B) = scale*C or
     !> op(A)*X - X*op(B) = scale*C,
     !> where op(A) = A or A**H, and A and B are both upper triangular. A is
     !> M-by-M and B is N-by-N; the right hand side C and the solution X are
     !> M-by-N; and scale is an output scale factor, set <= 1 to avoid
     !> overflow in X.

     subroutine la_wtrsyl(trana,tranb,isgn,m,n,a,lda,b,ldb,c,ldc,scale,info)
        use la_constants_qp,only:one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trana,tranb
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,lda,ldb,ldc,m,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notrna,notrnb
           integer(ilp) :: j,k,l
           real(qp) :: bignum,da11,db,eps,scaloc,sgn,smin,smlnum
           complex(qp) :: a11,suml,sumr,vec,x11
           ! Local Arrays
           real(qp) :: dum(1)
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Executable Statements
           ! decode and test input parameters
           notrna = la_lsame(trana,'N')
           notrnb = la_lsame(tranb,'N')
           info = 0
           if (.not. notrna .and. .not. la_lsame(trana,'C')) then
              info = -1
           else if (.not. notrnb .and. .not. la_lsame(tranb,'C')) then
              info = -2
           else if (isgn /= 1 .and. isgn /= -1) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,m)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WTRSYL',-info)
              return
           end if
           ! quick return if possible
           scale = one
           if (m == 0 .or. n == 0) return
           ! set constants to control overflow
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = smlnum*real(m*n,KIND=qp)/eps
           bignum = one/smlnum
           smin = max(smlnum,eps*la_wlange('M',m,m,a,lda,dum),eps*la_wlange('M', &
                      n,n,b,ldb,dum))
           sgn = isgn
           if (notrna .and. notrnb) then
              ! solve    a*x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                  ! a(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                        l-1
                ! r(k,l) = sum [a(k,i)*x(i,l)] +isgn*sum [x(k,j)*b(j,l)].
                        ! i=k+1                      j=1
              loop_30: do l = 1,n
                 do k = m,1,-1
                    suml = la_wdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_wdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = a(k,k) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=qp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=qp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_wladiv(vec*cmplx(scaloc,KIND=qp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_wqscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_30
           else if (.not. notrna .and. notrnb) then
              ! solve    a**h *x + isgn*x*b = scale*c.
              ! the (k,l)th block of x is determined starting from
              ! upper-left corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b(l,l) = c(k,l) - r(k,l)
              ! where
                         ! k-1                           l-1
                ! r(k,l) = sum [a**h(i,k)*x(i,l)] + isgn*sum [x(k,j)*b(j,l)]
                         ! i=1                           j=1
              loop_60: do l = 1,n
                 do k = 1,m
                    suml = la_wdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_wdotu(l - 1,c(k,1),ldc,b(1,l),1)
                    vec = c(k,l) - (suml + sgn*sumr)
                    scaloc = one
                    a11 = conjg(a(k,k)) + sgn*b(l,l)
                    da11 = abs(real(a11,KIND=qp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=qp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_wladiv(vec*cmplx(scaloc,KIND=qp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_wqscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_60
           else if (.not. notrna .and. .not. notrnb) then
              ! solve    a**h*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! upper-right corner column by column by
                  ! a**h(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! k-1
                 ! r(k,l) = sum [a**h(i,k)*x(i,l)] +
                          ! i=1
                                 ! n
                           ! isgn*sum [x(k,j)*b**h(l,j)].
                                ! j=l+1
              loop_90: do l = n,1,-1
                 do k = 1,m
                    suml = la_wdotc(k - 1,a(1,k),1,c(1,l),1)
                    sumr = la_wdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = conjg(a(k,k) + sgn*b(l,l))
                    da11 = abs(real(a11,KIND=qp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=qp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_wladiv(vec*cmplx(scaloc,KIND=qp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_wqscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_90
           else if (notrna .and. .not. notrnb) then
              ! solve    a*x + isgn*x*b**h = c.
              ! the (k,l)th block of x is determined starting from
              ! bottom-left corner column by column by
                 ! a(k,k)*x(k,l) + isgn*x(k,l)*b**h(l,l) = c(k,l) - r(k,l)
              ! where
                          ! m                          n
                ! r(k,l) = sum [a(k,i)*x(i,l)] + isgn*sum [x(k,j)*b**h(l,j)]
                        ! i=k+1                      j=l+1
              loop_120: do l = n,1,-1
                 do k = m,1,-1
                    suml = la_wdotu(m - k,a(k,min(k + 1,m)),lda,c(min(k + 1,m),l),1 &
                              )
                    sumr = la_wdotc(n - l,c(k,min(l + 1,n)),ldc,b(l,min(l + 1,n)), &
                              ldb)
                    vec = c(k,l) - (suml + sgn*conjg(sumr))
                    scaloc = one
                    a11 = a(k,k) + sgn*conjg(b(l,l))
                    da11 = abs(real(a11,KIND=qp)) + abs(aimag(a11))
                    if (da11 <= smin) then
                       a11 = smin
                       da11 = smin
                       info = 1
                    end if
                    db = abs(real(vec,KIND=qp)) + abs(aimag(vec))
                    if (da11 < one .and. db > one) then
                       if (db > bignum*da11) scaloc = one/db
                    end if
                    x11 = la_wladiv(vec*cmplx(scaloc,KIND=qp),a11)
                    if (scaloc /= one) then
                       do j = 1,n
                          call la_wqscal(m,scaloc,c(1,j),1)
                       end do
                       scale = scale*scaloc
                    end if
                    c(k,l) = x11
                 end do
              end do loop_120
           end if
           return
     end subroutine la_wtrsyl
#endif

     !> CHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a complex upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_chsein(side,eigsrc,initv,select,n,h,ldh,w,vl,ldvl,vr,ldvr,mm, &
               m,work,rwork,ifaill,ifailr,info)
        use la_constants_sp,only:czero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: h(ldh,*)
           complex(sp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),w(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: rzero = 0.0e+0_sp

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ks,ldwork
           real(sp) :: eps3,hnorm,smlnum,ulp,unfl
           complex(sp) :: cdum,wk
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -12
           else if (mm < m) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_slamch('SAFE MINIMUM')
           ulp = la_slamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ldwork = n
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ks = 1
           loop_100: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are czero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == czero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == czero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_clanhs('I',kr - kl + 1,h(kl,kl),ldh,rwork)
                    if (la_sisnan(hnorm)) then
                       info = -6
                       return
                    else if ((hnorm > rzero)) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wk = w(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. cabs1(w(i) - wk) < eps3) then
                       wk = wk + eps3
                       go to 60
                    end if
                 end do
                 w(k) = wk
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_claein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wk,vl(kl,ks) &
                              ,work,ldwork,rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifaill(ks) = k
                    else
                       ifaill(ks) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ks) = czero
                    end do
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_claein(.true.,noinit,kr,h,ldh,wk,vr(1,ks),work,ldwork, &
                              rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifailr(ks) = k
                    else
                       ifailr(ks) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ks) = czero
                    end do
                 end if
                 ks = ks + 1
              end if
           end do loop_100
           return
     end subroutine la_chsein
     !> ZHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a complex upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_zhsein(side,eigsrc,initv,select,n,h,ldh,w,vl,ldvl,vr,ldvr,mm, &
               m,work,rwork,ifaill,ifailr,info)
        use la_constants_dp,only:czero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: h(ldh,*)
           complex(dp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),w(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: rzero = 0.0e+0_dp

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ks,ldwork
           real(dp) :: eps3,hnorm,smlnum,ulp,unfl
           complex(dp) :: cdum,wk
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -12
           else if (mm < m) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_dlamch('SAFE MINIMUM')
           ulp = la_dlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ldwork = n
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ks = 1
           loop_100: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are czero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == czero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == czero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_zlanhs('I',kr - kl + 1,h(kl,kl),ldh,rwork)
                    if (la_disnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > rzero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wk = w(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. cabs1(w(i) - wk) < eps3) then
                       wk = wk + eps3
                       go to 60
                    end if
                 end do
                 w(k) = wk
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_zlaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wk,vl(kl,ks) &
                              ,work,ldwork,rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifaill(ks) = k
                    else
                       ifaill(ks) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ks) = czero
                    end do
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_zlaein(.true.,noinit,kr,h,ldh,wk,vr(1,ks),work,ldwork, &
                              rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifailr(ks) = k
                    else
                       ifailr(ks) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ks) = czero
                    end do
                 end if
                 ks = ks + 1
              end if
           end do loop_100
           return
     end subroutine la_zhsein
#ifdef LA_WITH_XDP
     !> YHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a complex upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_yhsein(side,eigsrc,initv,select,n,h,ldh,w,vl,ldvl,vr,ldvr,mm, &
               m,work,rwork,ifaill,ifailr,info)
        use la_constants_xdp,only:czero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(in) :: h(ldh,*)
           complex(xdp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),w(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: rzero = 0.0e+0_xdp

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ks,ldwork
           real(xdp) :: eps3,hnorm,smlnum,ulp,unfl
           complex(xdp) :: cdum,wk
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -12
           else if (mm < m) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('YHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_xlamch('SAFE MINIMUM')
           ulp = la_xlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ldwork = n
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ks = 1
           loop_100: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are czero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == czero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == czero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_ylanhs('I',kr - kl + 1,h(kl,kl),ldh,rwork)
                    if (la_xisnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > rzero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wk = w(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. cabs1(w(i) - wk) < eps3) then
                       wk = wk + eps3
                       go to 60
                    end if
                 end do
                 w(k) = wk
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_ylaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wk,vl(kl,ks) &
                              ,work,ldwork,rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifaill(ks) = k
                    else
                       ifaill(ks) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ks) = czero
                    end do
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_ylaein(.true.,noinit,kr,h,ldh,wk,vr(1,ks),work,ldwork, &
                              rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifailr(ks) = k
                    else
                       ifailr(ks) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ks) = czero
                    end do
                 end if
                 ks = ks + 1
              end if
           end do loop_100
           return
     end subroutine la_yhsein
#endif
#ifdef LA_WITH_QP
     !> WHSEIN: uses inverse iteration to find specified right and/or left
     !> eigenvectors of a complex upper Hessenberg matrix H.
     !> The right eigenvector x and the left eigenvector y of the matrix H
     !> corresponding to an eigenvalue w are defined by:
     !> H * x = w * x,     y**h * H = w * y**h
     !> where y**h denotes the conjugate transpose of the vector y.

     subroutine la_whsein(side,eigsrc,initv,select,n,h,ldh,w,vl,ldvl,vr,ldvr,mm, &
               m,work,rwork,ifaill,ifailr,info)
        use la_constants_qp,only:czero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: eigsrc,initv,side
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldh,ldvl,ldvr,mm,n
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           integer(ilp),intent(out) :: ifaill(*),ifailr(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: h(ldh,*)
           complex(qp),intent(inout) :: vl(ldvl,*),vr(ldvr,*),w(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: rzero = 0.0e+0_qp

           ! Local Scalars
           logical(lk) :: bothv,fromqr,leftv,noinit,rightv
           integer(ilp) :: i,iinfo,k,kl,kln,kr,ks,ldwork
           real(qp) :: eps3,hnorm,smlnum,ulp,unfl
           complex(qp) :: cdum,wk
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! decode and test the input parameters.
           bothv = la_lsame(side,'B')
           rightv = la_lsame(side,'R') .or. bothv
           leftv = la_lsame(side,'L') .or. bothv
           fromqr = la_lsame(eigsrc,'Q')
           noinit = la_lsame(initv,'N')
           ! set m to the number of columns required to store the selected
           ! eigenvectors.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           info = 0
           if (.not. rightv .and. .not. leftv) then
              info = -1
           else if (.not. fromqr .and. .not. la_lsame(eigsrc,'N')) then
              info = -2
           else if (.not. noinit .and. .not. la_lsame(initv,'U')) then
              info = -3
           else if (n < 0) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
              info = -10
           else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
              info = -12
           else if (mm < m) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WHSEIN',-info)
              return
           end if
           ! quick return if possible.
           if (n == 0) return
           ! set machine-dependent constants.
           unfl = la_qlamch('SAFE MINIMUM')
           ulp = la_qlamch('PRECISION')
           smlnum = unfl*(n/ulp)
           ldwork = n
           kl = 1
           kln = 0
           if (fromqr) then
              kr = 0
           else
              kr = n
           end if
           ks = 1
           loop_100: do k = 1,n
              if (select(k)) then
                 ! compute eigenvector(s) corresponding to w(k).
                 if (fromqr) then
                    ! if affiliation of eigenvalues is known, check whether
                    ! the matrix splits.
                    ! determine kl and kr such that 1 <= kl <= k <= kr <= n
                    ! and h(kl,kl-1) and h(kr+1,kr) are czero (or kl = 1 or
                    ! kr = n).
                    ! then inverse iteration can be performed with the
                    ! submatrix h(kl:n,kl:n) for a left eigenvector, and with
                    ! the submatrix h(1:kr,1:kr) for a right eigenvector.
                    do i = k,kl + 1,-1
                       if (h(i,i - 1) == czero) go to 30
                    end do
                    30 continue
                    kl = i
                    if (k > kr) then
                       do i = k,n - 1
                          if (h(i + 1,i) == czero) go to 50
                       end do
                       50 continue
                       kr = i
                    end if
                 end if
                 if (kl /= kln) then
                    kln = kl
                    ! compute infinity-norm of submatrix h(kl:kr,kl:kr) if it
                    ! has not ben computed before.
                    hnorm = la_wlanhs('I',kr - kl + 1,h(kl,kl),ldh,rwork)
                    if (la_qisnan(hnorm)) then
                       info = -6
                       return
                    else if (hnorm > rzero) then
                       eps3 = hnorm*ulp
                    else
                       eps3 = smlnum
                    end if
                 end if
                 ! perturb eigenvalue if it is close to any previous
                 ! selected eigenvalues affiliated to the submatrix
                 ! h(kl:kr,kl:kr). close roots are modified by eps3.
                 wk = w(k)
                 60 continue
                 do i = k - 1,kl,-1
                    if (select(i) .and. cabs1(w(i) - wk) < eps3) then
                       wk = wk + eps3
                       go to 60
                    end if
                 end do
                 w(k) = wk
                 if (leftv) then
                    ! compute left eigenvector.
                    call la_wlaein(.false.,noinit,n - kl + 1,h(kl,kl),ldh,wk,vl(kl,ks) &
                              ,work,ldwork,rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifaill(ks) = k
                    else
                       ifaill(ks) = 0
                    end if
                    do i = 1,kl - 1
                       vl(i,ks) = czero
                    end do
                 end if
                 if (rightv) then
                    ! compute right eigenvector.
                    call la_wlaein(.true.,noinit,kr,h,ldh,wk,vr(1,ks),work,ldwork, &
                              rwork,eps3,smlnum,iinfo)
                    if (iinfo > 0) then
                       info = info + 1
                       ifailr(ks) = k
                    else
                       ifailr(ks) = 0
                    end if
                    do i = kr + 1,n
                       vr(i,ks) = czero
                    end do
                 end if
                 ks = ks + 1
              end if
           end do loop_100
           return
     end subroutine la_whsein
#endif

     !> CTRSEN: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
     !> the leading positions on the diagonal of the upper triangular matrix
     !> T, and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.

     subroutine la_ctrsen(job,compq,select,n,t,ldt,q,ldq,w,m,s,sep,work,lwork, &
               info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,lwork,n
           real(sp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           complex(sp),intent(inout) :: q(ldq,*),t(ldt,*)
           complex(sp),intent(out) :: w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,ks,lwmin,n1,n2,nn
           real(sp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(sp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode and test the input parameters.
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           ! set m to the number of selected eigenvalues.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           n1 = m
           n2 = n - m
           nn = n1*n2
           info = 0
           lquery = (lwork == -1)
           if (wantsp) then
              lwmin = max(1,2*nn)
           else if (la_lsame(job,'N')) then
              lwmin = 1
           else if (la_lsame(job,'E')) then
              lwmin = max(1,nn)
           end if
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else if (lwork < lwmin .and. .not. lquery) then
              info = -14
           end if
           if (info == 0) then
              work(1) = lwmin
           end if
           if (info /= 0) then
              call la_xerbla('CTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_clange('1',n,n,t,ldt,rwork)
              go to 40
           end if
           ! collect the selected eigenvalues at the top left corner of t.
           ks = 0
           do k = 1,n
              if (select(k)) then
                 ks = ks + 1
                 ! swap the k-th eigenvalue to position ks.
                 if (k /= ks) call la_ctrexc(compq,n,t,ldt,q,ldq,k,ks,ierr)
              end if
           end do
           if (wants) then
              ! solve the sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_clacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_ctrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_clange('F',n1,n2,work,n1,rwork)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_clacn2(nn,work(nn + 1),work,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve t11*r - r*t22 = scale*x.
                    call la_ctrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**h*r - r*t22**h = scale*x.
                    call la_ctrsyl('C','C',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! copy reordered eigenvalues to w.
           do k = 1,n
              w(k) = t(k,k)
           end do
           work(1) = lwmin
           return
     end subroutine la_ctrsen
     !> ZTRSEN: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
     !> the leading positions on the diagonal of the upper triangular matrix
     !> T, and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.

     subroutine la_ztrsen(job,compq,select,n,t,ldt,q,ldq,w,m,s,sep,work,lwork, &
               info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,lwork,n
           real(dp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           complex(dp),intent(inout) :: q(ldq,*),t(ldt,*)
           complex(dp),intent(out) :: w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,ks,lwmin,n1,n2,nn
           real(dp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(dp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode and test the input parameters.
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           ! set m to the number of selected eigenvalues.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           n1 = m
           n2 = n - m
           nn = n1*n2
           info = 0
           lquery = (lwork == -1)
           if (wantsp) then
              lwmin = max(1,2*nn)
           else if (la_lsame(job,'N')) then
              lwmin = 1
           else if (la_lsame(job,'E')) then
              lwmin = max(1,nn)
           end if
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else if (lwork < lwmin .and. .not. lquery) then
              info = -14
           end if
           if (info == 0) then
              work(1) = lwmin
           end if
           if (info /= 0) then
              call la_xerbla('ZTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_zlange('1',n,n,t,ldt,rwork)
              go to 40
           end if
           ! collect the selected eigenvalues at the top left corner of t.
           ks = 0
           do k = 1,n
              if (select(k)) then
                 ks = ks + 1
                 ! swap the k-th eigenvalue to position ks.
                 if (k /= ks) call la_ztrexc(compq,n,t,ldt,q,ldq,k,ks,ierr)
              end if
           end do
           if (wants) then
              ! solve the sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_zlacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_ztrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_zlange('F',n1,n2,work,n1,rwork)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_zlacn2(nn,work(nn + 1),work,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve t11*r - r*t22 = scale*x.
                    call la_ztrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**h*r - r*t22**h = scale*x.
                    call la_ztrsyl('C','C',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! copy reordered eigenvalues to w.
           do k = 1,n
              w(k) = t(k,k)
           end do
           work(1) = lwmin
           return
     end subroutine la_ztrsen
#ifdef LA_WITH_XDP
     !> YTRSEN: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
     !> the leading positions on the diagonal of the upper triangular matrix
     !> T, and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.

     subroutine la_ytrsen(job,compq,select,n,t,ldt,q,ldq,w,m,s,sep,work,lwork, &
               info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,lwork,n
           real(xdp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           complex(xdp),intent(inout) :: q(ldq,*),t(ldt,*)
           complex(xdp),intent(out) :: w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,ks,lwmin,n1,n2,nn
           real(xdp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(xdp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode and test the input parameters.
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           ! set m to the number of selected eigenvalues.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           n1 = m
           n2 = n - m
           nn = n1*n2
           info = 0
           lquery = (lwork == -1)
           if (wantsp) then
              lwmin = max(1,2*nn)
           else if (la_lsame(job,'N')) then
              lwmin = 1
           else if (la_lsame(job,'E')) then
              lwmin = max(1,nn)
           end if
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else if (lwork < lwmin .and. .not. lquery) then
              info = -14
           end if
           if (info == 0) then
              work(1) = lwmin
           end if
           if (info /= 0) then
              call la_xerbla('YTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_ylange('1',n,n,t,ldt,rwork)
              go to 40
           end if
           ! collect the selected eigenvalues at the top left corner of t.
           ks = 0
           do k = 1,n
              if (select(k)) then
                 ks = ks + 1
                 ! swap the k-th eigenvalue to position ks.
                 if (k /= ks) call la_ytrexc(compq,n,t,ldt,q,ldq,k,ks,ierr)
              end if
           end do
           if (wants) then
              ! solve the sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_ylacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_ytrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_ylange('F',n1,n2,work,n1,rwork)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_ylacn2(nn,work(nn + 1),work,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve t11*r - r*t22 = scale*x.
                    call la_ytrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**h*r - r*t22**h = scale*x.
                    call la_ytrsyl('C','C',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! copy reordered eigenvalues to w.
           do k = 1,n
              w(k) = t(k,k)
           end do
           work(1) = lwmin
           return
     end subroutine la_ytrsen
#endif
#ifdef LA_WITH_QP
     !> WTRSEN: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that a selected cluster of eigenvalues appears in
     !> the leading positions on the diagonal of the upper triangular matrix
     !> T, and the leading columns of Q form an orthonormal basis of the
     !> corresponding right invariant subspace.
     !> Optionally the routine computes the reciprocal condition numbers of
     !> the cluster of eigenvalues and/or the invariant subspace.

     subroutine la_wtrsen(job,compq,select,n,t,ldt,q,ldq,w,m,s,sep,work,lwork, &
               info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,job
           integer(ilp),intent(out) :: info,m
           integer(ilp),intent(in) :: ldq,ldt,lwork,n
           real(qp),intent(out) :: s,sep
           ! Array Arguments
           logical(lk),intent(in) :: select(*)
           complex(qp),intent(inout) :: q(ldq,*),t(ldt,*)
           complex(qp),intent(out) :: w(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantbh,wantq,wants,wantsp
           integer(ilp) :: ierr,k,kase,ks,lwmin,n1,n2,nn
           real(qp) :: est,rnorm,scale
           ! Local Arrays
           integer(ilp) :: isave(3)
           real(qp) :: rwork(1)
           ! Intrinsic Functions
           intrinsic :: max,sqrt
           ! Executable Statements
           ! decode and test the input parameters.
           wantbh = la_lsame(job,'B')
           wants = la_lsame(job,'E') .or. wantbh
           wantsp = la_lsame(job,'V') .or. wantbh
           wantq = la_lsame(compq,'V')
           ! set m to the number of selected eigenvalues.
           m = 0
           do k = 1,n
              if (select(k)) m = m + 1
           end do
           n1 = m
           n2 = n - m
           nn = n1*n2
           info = 0
           lquery = (lwork == -1)
           if (wantsp) then
              lwmin = max(1,2*nn)
           else if (la_lsame(job,'N')) then
              lwmin = 1
           else if (la_lsame(job,'E')) then
              lwmin = max(1,nn)
           end if
           if (.not. la_lsame(job,'N') .and. .not. wants .and. .not. wantsp) then
              info = -1
           else if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -2
           else if (n < 0) then
              info = -4
           else if (ldt < max(1,n)) then
              info = -6
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -8
           else if (lwork < lwmin .and. .not. lquery) then
              info = -14
           end if
           if (info == 0) then
              work(1) = lwmin
           end if
           if (info /= 0) then
              call la_xerbla('WTRSEN',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == n .or. m == 0) then
              if (wants) s = one
              if (wantsp) sep = la_wlange('1',n,n,t,ldt,rwork)
              go to 40
           end if
           ! collect the selected eigenvalues at the top left corner of t.
           ks = 0
           do k = 1,n
              if (select(k)) then
                 ks = ks + 1
                 ! swap the k-th eigenvalue to position ks.
                 if (k /= ks) call la_wtrexc(compq,n,t,ldt,q,ldq,k,ks,ierr)
              end if
           end do
           if (wants) then
              ! solve the sylvester equation for r:
                 ! t11*r - r*t22 = scale*t12
              call la_wlacpy('F',n1,n2,t(1,n1 + 1),ldt,work,n1)
              call la_wtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work,n1, &
                        scale,ierr)
              ! estimate the reciprocal of the condition number of the cluster
              ! of eigenvalues.
              rnorm = la_wlange('F',n1,n2,work,n1,rwork)
              if (rnorm == zero) then
                 s = one
              else
                 s = scale/(sqrt(scale*scale/rnorm + rnorm)*sqrt(rnorm))
              end if
           end if
           if (wantsp) then
              ! estimate sep(t11,t22).
              est = zero
              kase = 0
              30 continue
              call la_wlacn2(nn,work(nn + 1),work,est,kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! solve t11*r - r*t22 = scale*x.
                    call la_wtrsyl('N','N',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 else
                    ! solve t11**h*r - r*t22**h = scale*x.
                    call la_wtrsyl('C','C',-1,n1,n2,t,ldt,t(n1 + 1,n1 + 1),ldt,work, &
                              n1,scale,ierr)
                 end if
                 go to 30
              end if
              sep = scale/est
           end if
           40 continue
           ! copy reordered eigenvalues to w.
           do k = 1,n
              w(k) = t(k,k)
           end do
           work(1) = lwmin
           return
     end subroutine la_wtrsen
#endif

     !> CHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**H, where T is an upper triangular matrix (the
     !> Schur form), and Z is the unitary matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input unitary
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.

     pure subroutine la_chseqr(job,compz,n,ilo,ihi,h,ldh,w,z,ldz,work,lwork,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           complex(sp),intent(inout) :: h(ldh,*),z(ldz,*)
           complex(sp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           real(sp),parameter :: rzero = 0.0_sp
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_clahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_clahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           complex(sp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: cmplx,max,min,real
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = cmplx(real(max(1,n),KIND=sp),rzero,KIND=sp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -10
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('CHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_claqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                        lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(work(1),KIND=sp),real(max(1,n),KIND=sp)), &
                        rzero,KIND=sp)
              return
           else
              ! ==== copy eigenvalues isolated by la_cgebal ====
              if (ilo > 1) call la_ccopy(ilo - 1,h,ldh + 1,w,1)
              if (ihi < n) call la_ccopy(n - ihi,h(ihi + 1,ihi + 1),ldh + 1,w(ihi + 1),1)

              ! ==== initialize z, if requested ====
              if (initz) call la_claset('A',n,n,czero,cone,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 w(ilo) = h(ilo,ilo)
                 return
              end if
              ! ==== la_clahqr/la_claqr0 crossover point ====
              nmin = la_ilaenv(12,'CHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_claqr0 for big matrices; la_clahqr for small ones ====
              if (n > nmin) then
                 call la_claqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                           lwork,info)
              else
                 ! ==== small matrix ====
                 call la_clahqr(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,info)

                 if (info > 0) then
                    ! ==== a rare la_clahqr failure!  la_claqr0 sometimes succeeds
                    ! .    when la_clahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_claqr0 directly. ====
                       call la_claqr0(wantt,wantz,n,ilo,kbot,h,ldh,w,ilo,ihi,z,ldz, &
                                  work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_claqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_claqr0. ====
                       call la_clacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = czero
                       call la_claset('A',nl,nl - n,czero,czero,hl(1,n + 1),nl)
                       call la_claqr0(wantt,wantz,nl,ilo,kbot,hl,nl,w,ilo,ihi,z, &
                                 ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_clacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_claset('L',n - 2,n - 2,czero, &
                        czero,h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(max(1,n),KIND=sp),real(work(1),KIND=sp)), &
                        rzero,KIND=sp)
           end if
     end subroutine la_chseqr
     !> ZHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**H, where T is an upper triangular matrix (the
     !> Schur form), and Z is the unitary matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input unitary
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.

     pure subroutine la_zhseqr(job,compz,n,ilo,ihi,h,ldh,w,z,ldz,work,lwork,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           complex(dp),intent(inout) :: h(ldh,*),z(ldz,*)
           complex(dp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           real(dp),parameter :: rzero = 0.0_dp
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_zlahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_zlahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           complex(dp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = cmplx(real(max(1,n),KIND=dp),rzero,KIND=dp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -10
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('ZHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_zlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                        lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(work(1),KIND=dp),real(max(1,n),KIND=dp)), &
                        rzero,KIND=dp)
              return
           else
              ! ==== copy eigenvalues isolated by la_zgebal ====
              if (ilo > 1) call la_zcopy(ilo - 1,h,ldh + 1,w,1)
              if (ihi < n) call la_zcopy(n - ihi,h(ihi + 1,ihi + 1),ldh + 1,w(ihi + 1),1)

              ! ==== initialize z, if requested ====
              if (initz) call la_zlaset('A',n,n,czero,cone,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 w(ilo) = h(ilo,ilo)
                 return
              end if
              ! ==== la_zlahqr/la_zlaqr0 crossover point ====
              nmin = la_ilaenv(12,'ZHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_zlaqr0 for big matrices; la_zlahqr for small ones ====
              if (n > nmin) then
                 call la_zlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                           lwork,info)
              else
                 ! ==== small matrix ====
                 call la_zlahqr(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,info)

                 if (info > 0) then
                    ! ==== a rare la_zlahqr failure!  la_zlaqr0 sometimes succeeds
                    ! .    when la_zlahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_zlaqr0 directly. ====
                       call la_zlaqr0(wantt,wantz,n,ilo,kbot,h,ldh,w,ilo,ihi,z,ldz, &
                                  work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_zlaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_zlaqr0. ====
                       call la_zlacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = czero
                       call la_zlaset('A',nl,nl - n,czero,czero,hl(1,n + 1),nl)
                       call la_zlaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,w,ilo,ihi,z, &
                                 ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_zlacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_zlaset('L',n - 2,n - 2,czero, &
                        czero,h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(max(1,n),KIND=dp),real(work(1),KIND=dp)), &
                        rzero,KIND=dp)
           end if
     end subroutine la_zhseqr
#ifdef LA_WITH_XDP
     !> YHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**H, where T is an upper triangular matrix (the
     !> Schur form), and Z is the unitary matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input unitary
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.

     pure subroutine la_yhseqr(job,compz,n,ilo,ihi,h,ldh,w,z,ldz,work,lwork,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           complex(xdp),intent(inout) :: h(ldh,*),z(ldz,*)
           complex(xdp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           real(xdp),parameter :: rzero = 0.0_xdp
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_ylahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_ylahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           complex(xdp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = cmplx(real(max(1,n),KIND=xdp),rzero,KIND=xdp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -10
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('YHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_ylaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                        lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(work(1),KIND=xdp),real(max(1,n),KIND=xdp)), &
                        rzero,KIND=xdp)
              return
           else
              ! ==== copy eigenvalues isolated by la_ygebal ====
              if (ilo > 1) call la_ycopy(ilo - 1,h,ldh + 1,w,1)
              if (ihi < n) call la_ycopy(n - ihi,h(ihi + 1,ihi + 1),ldh + 1,w(ihi + 1),1)

              ! ==== initialize z, if requested ====
              if (initz) call la_ylaset('A',n,n,czero,cone,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 w(ilo) = h(ilo,ilo)
                 return
              end if
              ! ==== la_ylahqr/la_ylaqr0 crossover point ====
              nmin = la_ilaenv(12,'YHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_ylaqr0 for big matrices; la_ylahqr for small ones ====
              if (n > nmin) then
                 call la_ylaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                           lwork,info)
              else
                 ! ==== small matrix ====
                 call la_ylahqr(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,info)

                 if (info > 0) then
                    ! ==== a rare la_ylahqr failure!  la_ylaqr0 sometimes succeeds
                    ! .    when la_ylahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_ylaqr0 directly. ====
                       call la_ylaqr0(wantt,wantz,n,ilo,kbot,h,ldh,w,ilo,ihi,z,ldz, &
                                  work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_ylaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_ylaqr0. ====
                       call la_ylacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = czero
                       call la_ylaset('A',nl,nl - n,czero,czero,hl(1,n + 1),nl)
                       call la_ylaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,w,ilo,ihi,z, &
                                 ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_ylacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_ylaset('L',n - 2,n - 2,czero, &
                        czero,h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(max(1,n),KIND=xdp),real(work(1),KIND=xdp)), &
                        rzero,KIND=xdp)
           end if
     end subroutine la_yhseqr
#endif
#ifdef LA_WITH_QP
     !> WHSEQR: computes the eigenvalues of a Hessenberg matrix H
     !> and, optionally, the matrices T and Z from the Schur decomposition
     !> H = Z T Z**H, where T is an upper triangular matrix (the
     !> Schur form), and Z is the unitary matrix of Schur vectors.
     !> Optionally Z may be postmultiplied into an input unitary
     !> matrix Q so that this routine can give the Schur factorization
     !> of a matrix A which has been reduced to the Hessenberg form H
     !> by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.

     pure subroutine la_whseqr(job,compz,n,ilo,ihi,h,ldh,w,z,ldz,work,lwork,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ldh,ldz,lwork,n
           integer(ilp),intent(out) :: info
           character,intent(in) :: compz,job
           ! Array Arguments
           complex(qp),intent(inout) :: h(ldh,*),z(ldz,*)
           complex(qp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ntiny = 15
           integer(ilp),parameter :: nl = 49
           real(qp),parameter :: rzero = 0.0_qp
           ! ==== matrices of order ntiny or smaller must be processed by
           ! .    la_wlahqr because of insufficient subdiagonal scratch space.
           ! .    (this is a hard limit.) ====

           ! ==== nl allocates some local workspace to help small matrices
           ! .    through a rare la_wlahqr failure.  nl > ntiny = 15 is
           ! .    required and nl <= nmin = la_ilaenv(ispec=12,...) is recom-
           ! .    mended.  (the default value of nmin is 75.)  using nl = 49
           ! .    allows up to six simultaneous shifts and a 16-by-16
           ! .    deflation window.  ====

           ! Local Arrays
           complex(qp) :: hl(nl,nl),workl(nl)
           ! Local Scalars
           integer(ilp) :: kbot,nmin
           logical(lk) :: initz,lquery,wantt,wantz
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max,min
           ! Executable Statements
           ! ==== decode and check the input parameters. ====
           wantt = la_lsame(job,'S')
           initz = la_lsame(compz,'I')
           wantz = initz .or. la_lsame(compz,'V')
           work(1) = cmplx(real(max(1,n),KIND=qp),rzero,KIND=qp)
           lquery = lwork == -1
           info = 0
           if (.not. la_lsame(job,'E') .and. .not. wantt) then
              info = -1
           else if (.not. la_lsame(compz,'N') .and. .not. wantz) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ilo < 1 .or. ilo > max(1,n)) then
              info = -4
           else if (ihi < min(ilo,n) .or. ihi > n) then
              info = -5
           else if (ldh < max(1,n)) then
              info = -7
           else if (ldz < 1 .or. (wantz .and. ldz < max(1,n))) then
              info = -10
           else if (lwork < max(1,n) .and. .not. lquery) then
              info = -12
           end if
           if (info /= 0) then
              ! ==== quick return in case of invalid argument. ====
              call la_xerbla('WHSEQR',-info)
              return
           else if (n == 0) then
              ! ==== quick return in case n = 0; nothing to do. ====
              return
           else if (lquery) then
              ! ==== quick return in case of a workspace query ====
              call la_wlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                        lwork,info)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(work(1),KIND=qp),real(max(1,n),KIND=qp)), &
                        rzero,KIND=qp)
              return
           else
              ! ==== copy eigenvalues isolated by la_wgebal ====
              if (ilo > 1) call la_wcopy(ilo - 1,h,ldh + 1,w,1)
              if (ihi < n) call la_wcopy(n - ihi,h(ihi + 1,ihi + 1),ldh + 1,w(ihi + 1),1)

              ! ==== initialize z, if requested ====
              if (initz) call la_wlaset('A',n,n,czero,cone,z,ldz)
              ! ==== quick return if possible ====
              if (ilo == ihi) then
                 w(ilo) = h(ilo,ilo)
                 return
              end if
              ! ==== la_wlahqr/la_wlaqr0 crossover point ====
              nmin = la_ilaenv(12,'WHSEQR',job(:1)//compz(:1),n,ilo,ihi,lwork)

              nmin = max(ntiny,nmin)
              ! ==== la_wlaqr0 for big matrices; la_wlahqr for small ones ====
              if (n > nmin) then
                 call la_wlaqr0(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,work, &
                           lwork,info)
              else
                 ! ==== small matrix ====
                 call la_wlahqr(wantt,wantz,n,ilo,ihi,h,ldh,w,ilo,ihi,z,ldz,info)

                 if (info > 0) then
                    ! ==== a rare la_wlahqr failure!  la_wlaqr0 sometimes succeeds
                    ! .    when la_wlahqr fails. ====
                    kbot = info
                    if (n >= nl) then
                       ! ==== larger matrices have enough subdiagonal scratch
                       ! .    space to call la_wlaqr0 directly. ====
                       call la_wlaqr0(wantt,wantz,n,ilo,kbot,h,ldh,w,ilo,ihi,z,ldz, &
                                  work,lwork,info)
                    else
                       ! ==== tiny matrices don't have enough subdiagonal
                       ! .    scratch space to benefit from la_wlaqr0.  hence,
                       ! .    tiny matrices must be copied into a larger
                       ! .    array before calling la_wlaqr0. ====
                       call la_wlacpy('A',n,n,h,ldh,hl,nl)
                       hl(n + 1,n) = czero
                       call la_wlaset('A',nl,nl - n,czero,czero,hl(1,n + 1),nl)
                       call la_wlaqr0(wantt,wantz,nl,ilo,kbot,hl,nl,w,ilo,ihi,z, &
                                 ldz,workl,nl,info)
                       if (wantt .or. info /= 0) call la_wlacpy('A',n,n,hl,nl,h,ldh)

                    end if
                 end if
              end if
              ! ==== clear out the trash, if necessary. ====
              if ((wantt .or. info /= 0) .and. n > 2) call la_wlaset('L',n - 2,n - 2,czero, &
                        czero,h(3,1),ldh)
              ! ==== ensure reported workspace size is backward-compatible with
              ! .    previous lapack versions. ====
              work(1) = cmplx(max(real(max(1,n),KIND=qp),real(work(1),KIND=qp)), &
                        rzero,KIND=qp)
           end if
     end subroutine la_whseqr
#endif

end module la_lapack_eigv_gen2
