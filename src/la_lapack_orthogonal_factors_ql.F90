!> LQ and QL factorizations: blocked, short-wide and triangular-pentagonal variants
module la_lapack_orthogonal_factors_ql
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_householder_reflectors
     use la_lapack_orthogonal_factors_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sorg2l
     public :: la_sorgl2
     public :: la_sorglq
     public :: la_sorgql
     public :: la_sorm2l
     public :: la_sorml2
     public :: la_sormlq
     public :: la_sormql
     public :: la_sgemlqt
     public :: la_stplqt2
     public :: la_stpmlqt
     public :: la_sgelq2
     public :: la_sgelqf
     public :: la_sgelqt3
     public :: la_sgeql2
     public :: la_sgeqlf
     public :: la_slamswlq
     public :: la_stplqt
     public :: la_sgelqt
     public :: la_sgemlq
     public :: la_slaswlq
     public :: la_sgelq
     public :: la_dorg2l
     public :: la_dorgl2
     public :: la_dorglq
     public :: la_dorgql
     public :: la_dorm2l
     public :: la_dorml2
     public :: la_dormlq
     public :: la_dormql
     public :: la_dgemlqt
     public :: la_dtplqt2
     public :: la_dtpmlqt
     public :: la_dgelq2
     public :: la_dgelqf
     public :: la_dgelqt3
     public :: la_dgeql2
     public :: la_dgeqlf
     public :: la_dlamswlq
     public :: la_dtplqt
     public :: la_dgelqt
     public :: la_dgemlq
     public :: la_dlaswlq
     public :: la_dgelq
#ifdef LA_WITH_XDP
     public :: la_xorg2l
     public :: la_xorgl2
     public :: la_xorglq
     public :: la_xorgql
     public :: la_xorm2l
     public :: la_xorml2
     public :: la_xormlq
     public :: la_xormql
     public :: la_xgemlqt
     public :: la_xtplqt2
     public :: la_xtpmlqt
     public :: la_xgelq2
     public :: la_xgelqf
     public :: la_xgelqt3
     public :: la_xgeql2
     public :: la_xgeqlf
     public :: la_xlamswlq
     public :: la_xtplqt
     public :: la_xgelqt
     public :: la_xgemlq
     public :: la_xlaswlq
     public :: la_xgelq
#endif
#ifdef LA_WITH_QP
     public :: la_qorg2l
     public :: la_qorgl2
     public :: la_qorglq
     public :: la_qorgql
     public :: la_qorm2l
     public :: la_qorml2
     public :: la_qormlq
     public :: la_qormql
     public :: la_qgemlqt
     public :: la_qtplqt2
     public :: la_qtpmlqt
     public :: la_qgelq2
     public :: la_qgelqf
     public :: la_qgelqt3
     public :: la_qgeql2
     public :: la_qgeqlf
     public :: la_qlamswlq
     public :: la_qtplqt
     public :: la_qgelqt
     public :: la_qgemlq
     public :: la_qlaswlq
     public :: la_qgelq
#endif
     public :: la_ctplqt2
     public :: la_cung2l
     public :: la_cungl2
     public :: la_cunglq
     public :: la_cungql
     public :: la_cunm22
     public :: la_cunm2l
     public :: la_cunml2
     public :: la_cunmlq
     public :: la_cunmql
     public :: la_cgelq2
     public :: la_cgelqf
     public :: la_cgelqt3
     public :: la_cgemlqt
     public :: la_cgeql2
     public :: la_cgeqlf
     public :: la_ctplqt
     public :: la_ctpmlqt
     public :: la_cgelqt
     public :: la_clamswlq
     public :: la_claswlq
     public :: la_cgelq
     public :: la_cgemlq
     public :: la_ztplqt2
     public :: la_zung2l
     public :: la_zungl2
     public :: la_zunglq
     public :: la_zungql
     public :: la_zunm22
     public :: la_zunm2l
     public :: la_zunml2
     public :: la_zunmlq
     public :: la_zunmql
     public :: la_zgelq2
     public :: la_zgelqf
     public :: la_zgelqt3
     public :: la_zgemlqt
     public :: la_zgeql2
     public :: la_zgeqlf
     public :: la_ztplqt
     public :: la_ztpmlqt
     public :: la_zgelqt
     public :: la_zlamswlq
     public :: la_zlaswlq
     public :: la_zgelq
     public :: la_zgemlq
#ifdef LA_WITH_XDP
     public :: la_ytplqt2
     public :: la_yung2l
     public :: la_yungl2
     public :: la_yunglq
     public :: la_yungql
     public :: la_yunm22
     public :: la_yunm2l
     public :: la_yunml2
     public :: la_yunmlq
     public :: la_yunmql
     public :: la_ygelq2
     public :: la_ygelqf
     public :: la_ygelqt3
     public :: la_ygemlqt
     public :: la_ygeql2
     public :: la_ygeqlf
     public :: la_ytplqt
     public :: la_ytpmlqt
     public :: la_ygelqt
     public :: la_ylamswlq
     public :: la_ylaswlq
     public :: la_ygelq
     public :: la_ygemlq
#endif
#ifdef LA_WITH_QP
     public :: la_wtplqt2
     public :: la_wung2l
     public :: la_wungl2
     public :: la_wunglq
     public :: la_wungql
     public :: la_wunm22
     public :: la_wunm2l
     public :: la_wunml2
     public :: la_wunmlq
     public :: la_wunmql
     public :: la_wgelq2
     public :: la_wgelqf
     public :: la_wgelqt3
     public :: la_wgemlqt
     public :: la_wgeql2
     public :: la_wgeqlf
     public :: la_wtplqt
     public :: la_wtpmlqt
     public :: la_wgelqt
     public :: la_wlamswlq
     public :: la_wlaswlq
     public :: la_wgelq
     public :: la_wgemlq
#endif

     contains

     !> SORG2L: generates an m by n real matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by SGEQLF.

     pure subroutine la_sorg2l(m,n,k,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SORG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = zero
              end do
              a(m - n + j,j) = one
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = one
              call la_slarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_sscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = one - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to zero
              do l = m - n + ii + 1,m
                 a(l,ii) = zero
              end do
           end do
           return
     end subroutine la_sorg2l
     !> DORG2L: generates an m by n real matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by DGEQLF.

     pure subroutine la_dorg2l(m,n,k,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DORG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = zero
              end do
              a(m - n + j,j) = one
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = one
              call la_dlarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_dscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = one - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to zero
              do l = m - n + ii + 1,m
                 a(l,ii) = zero
              end do
           end do
           return
     end subroutine la_dorg2l
#ifdef LA_WITH_XDP
     !> XORG2L: generates an m by n real matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by XGEQLF.

     pure subroutine la_xorg2l(m,n,k,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('XORG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = zero
              end do
              a(m - n + j,j) = one
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = one
              call la_xlarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_xscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = one - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to zero
              do l = m - n + ii + 1,m
                 a(l,ii) = zero
              end do
           end do
           return
     end subroutine la_xorg2l
#endif
#ifdef LA_WITH_QP
     !> QORG2L: generates an m by n real matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by QGEQLF.

     pure subroutine la_qorg2l(m,n,k,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QORG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = zero
              end do
              a(m - n + j,j) = one
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = one
              call la_qlarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_qscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = one - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to zero
              do l = m - n + ii + 1,m
                 a(l,ii) = zero
              end do
           end do
           return
     end subroutine la_qorg2l
#endif

     !> SORGL2: generates an m by n real matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by SGELQF.

     pure subroutine la_sorgl2(m,n,k,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SORGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = zero
                 end do
                 if (j > k .and. j <= m) a(j,j) = one
              end do
           end if
           do i = k,1,-1
              ! apply h(i) to a(i:m,i:n) from the right
              if (i < n) then
                 if (i < m) then
                    a(i,i) = one
                    call la_slarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                              lda,work)
                 end if
                 call la_sscal(n - i,-tau(i),a(i,i + 1),lda)
              end if
              a(i,i) = one - tau(i)
              ! set a(i,1:i-1) to zero
              do l = 1,i - 1
                 a(i,l) = zero
              end do
           end do
           return
     end subroutine la_sorgl2
     !> DORGL2: generates an m by n real matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by DGELQF.

     pure subroutine la_dorgl2(m,n,k,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DORGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = zero
                 end do
                 if (j > k .and. j <= m) a(j,j) = one
              end do
           end if
           do i = k,1,-1
              ! apply h(i) to a(i:m,i:n) from the right
              if (i < n) then
                 if (i < m) then
                    a(i,i) = one
                    call la_dlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                              lda,work)
                 end if
                 call la_dscal(n - i,-tau(i),a(i,i + 1),lda)
              end if
              a(i,i) = one - tau(i)
              ! set a(i,1:i-1) to zero
              do l = 1,i - 1
                 a(i,l) = zero
              end do
           end do
           return
     end subroutine la_dorgl2
#ifdef LA_WITH_XDP
     !> XORGL2: generates an m by n real matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by XGELQF.

     pure subroutine la_xorgl2(m,n,k,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('XORGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = zero
                 end do
                 if (j > k .and. j <= m) a(j,j) = one
              end do
           end if
           do i = k,1,-1
              ! apply h(i) to a(i:m,i:n) from the right
              if (i < n) then
                 if (i < m) then
                    a(i,i) = one
                    call la_xlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                              lda,work)
                 end if
                 call la_xscal(n - i,-tau(i),a(i,i + 1),lda)
              end if
              a(i,i) = one - tau(i)
              ! set a(i,1:i-1) to zero
              do l = 1,i - 1
                 a(i,l) = zero
              end do
           end do
           return
     end subroutine la_xorgl2
#endif
#ifdef LA_WITH_QP
     !> QORGL2: generates an m by n real matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by QGELQF.

     pure subroutine la_qorgl2(m,n,k,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QORGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = zero
                 end do
                 if (j > k .and. j <= m) a(j,j) = one
              end do
           end if
           do i = k,1,-1
              ! apply h(i) to a(i:m,i:n) from the right
              if (i < n) then
                 if (i < m) then
                    a(i,i) = one
                    call la_qlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                              lda,work)
                 end if
                 call la_qscal(n - i,-tau(i),a(i,i + 1),lda)
              end if
              a(i,i) = one - tau(i)
              ! set a(i,1:i-1) to zero
              do l = 1,i - 1
                 a(i,l) = zero
              end do
           end do
           return
     end subroutine la_qorgl2
#endif

     !> SORGLQ: generates an M-by-N real matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by SGELQF.

     pure subroutine la_sorglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'SORGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('SORGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'SORGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'SORGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to zero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_sorgl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_slarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**t to a(i+ib:m,i:n) from the right
                    call la_slarfb('RIGHT','TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - i + &
                    1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
                 ! apply h**t to columns i:n of current block
                 call la_sorgl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to zero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_sorglq
     !> DORGLQ: generates an M-by-N real matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by DGELQF.

     pure subroutine la_dorglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'DORGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DORGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'DORGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'DORGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to zero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_dorgl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_dlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**t to a(i+ib:m,i:n) from the right
                    call la_dlarfb('RIGHT','TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - i + &
                    1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
                 ! apply h**t to columns i:n of current block
                 call la_dorgl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to zero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_dorglq
#ifdef LA_WITH_XDP
     !> XORGLQ: generates an M-by-N real matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by XGELQF.

     pure subroutine la_xorglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'XORGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('XORGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'XORGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'XORGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to zero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_xorgl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_xlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**t to a(i+ib:m,i:n) from the right
                    call la_xlarfb('RIGHT','TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - i + &
                    1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
                 ! apply h**t to columns i:n of current block
                 call la_xorgl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to zero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_xorglq
#endif
#ifdef LA_WITH_QP
     !> QORGLQ: generates an M-by-N real matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by QGELQF.

     pure subroutine la_qorglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'QORGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QORGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'QORGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'QORGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to zero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_qorgl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_qlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**t to a(i+ib:m,i:n) from the right
                    call la_qlarfb('RIGHT','TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - i + &
                    1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
                 ! apply h**t to columns i:n of current block
                 call la_qorgl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to zero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_qorglq
#endif

     !> SORGQL: generates an M-by-N real matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by SGEQLF.

     pure subroutine la_sorgql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'SORGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SORGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'SORGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'SORGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to zero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_sorg2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_slarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_slarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_sorg2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to zero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_sorgql
     !> DORGQL: generates an M-by-N real matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by DGEQLF.

     pure subroutine la_dorgql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'DORGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DORGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'DORGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'DORGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to zero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_dorg2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_dlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_dlarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_dorg2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to zero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_dorgql
#ifdef LA_WITH_XDP
     !> XORGQL: generates an M-by-N real matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by XGEQLF.

     pure subroutine la_xorgql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'XORGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XORGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'XORGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'XORGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to zero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_xorg2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_xlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_xlarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_xorg2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to zero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_xorgql
#endif
#ifdef LA_WITH_QP
     !> QORGQL: generates an M-by-N real matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by QGEQLF.

     pure subroutine la_qorgql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'QORGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QORGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'QORGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'QORGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to zero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = zero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_qorg2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_qlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_qlarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_qorg2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to zero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = zero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_qorgql
#endif

     !> SORM2L: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T * C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by SGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_sorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('SORM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i)
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = one
              call la_slarf(side,mi,ni,a(1,i),1,tau(i),c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_sorm2l
     !> DORM2L: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T * C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_dorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('DORM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i)
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = one
              call la_dlarf(side,mi,ni,a(1,i),1,tau(i),c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_dorm2l
#ifdef LA_WITH_XDP
     !> XORM2L: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T * C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by XGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_xorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           real(xdp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('XORM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i)
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = one
              call la_xlarf(side,mi,ni,a(1,i),1,tau(i),c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_xorm2l
#endif
#ifdef LA_WITH_QP
     !> QORM2L: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T * C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by QGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_qorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('QORM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i)
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = one
              call la_qlarf(side,mi,ni,a(1,i),1,tau(i),c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_qorm2l
#endif

     !> SORML2: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by SGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_sorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('SORML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i)
              aii = a(i,i)
              a(i,i) = one
              call la_slarf(side,mi,ni,a(i,i),lda,tau(i),c(ic,jc),ldc,work)

              a(i,i) = aii
           end do
           return
     end subroutine la_sorml2
     !> DORML2: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_dorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('DORML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i)
              aii = a(i,i)
              a(i,i) = one
              call la_dlarf(side,mi,ni,a(i,i),lda,tau(i),c(ic,jc),ldc,work)

              a(i,i) = aii
           end do
           return
     end subroutine la_dorml2
#ifdef LA_WITH_XDP
     !> XORML2: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by XGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_xorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           real(xdp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('XORML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i)
              aii = a(i,i)
              a(i,i) = one
              call la_xlarf(side,mi,ni,a(i,i),lda,tau(i),c(ic,jc),ldc,work)

              a(i,i) = aii
           end do
           return
     end subroutine la_xorml2
#endif
#ifdef LA_WITH_QP
     !> QORML2: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'T', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'T',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by QGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_qorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('QORML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i)
              aii = a(i,i)
              a(i,i) = one
              call la_qlarf(side,mi,ni,a(i,i),lda,tau(i),c(ic,jc),ldc,work)

              a(i,i) = aii
           end do
           return
     end subroutine la_qorml2
#endif

     !> SORMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_sormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'SORMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'SORMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_sorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_slarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**t is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**t
                 call la_slarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_sormlq
     !> DORMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_dormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'DORMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'DORMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_dorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_dlarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**t is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**t
                 call la_dlarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_dormlq
#ifdef LA_WITH_XDP
     !> XORMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by XGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_xormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'XORMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XORMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'XORMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_xorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_xlarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**t is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**t
                 call la_xlarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_xormlq
#endif
#ifdef LA_WITH_QP
     !> QORMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by QGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_qormlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'QORMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'QORMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_qorml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_qlarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**t is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**t
                 call la_qlarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_qormlq
#endif

     !> SORMQL: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by SGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_sormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'SORMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'SORMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_sorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_slarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**t is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**t
                 call la_slarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_sormql
     !> DORMQL: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_dormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'DORMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'DORMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_dorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_dlarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**t is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**t
                 call la_dlarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_dormql
#ifdef LA_WITH_XDP
     !> XORMQL: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by XGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_xormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'XORMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XORMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'XORMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_xorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_xlarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**t is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**t
                 call la_xlarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_xormql
#endif
#ifdef LA_WITH_QP
     !> QORMQL: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by QGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_qormql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'QORMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'QORMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_qorm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_qlarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**t is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**t is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**t
                 call la_qlarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_qormql
#endif

     !> DGEMLQT overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'T':   Q**T C            C Q**T
     !> where Q is a real orthogonal matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**T
     !> generated using the compact WY representation as returned by SGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_sgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           real(sp),intent(in) :: v(ldv,*),t(ldt,*)
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('SGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_slarfb('L','T','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_slarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_slarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_slarfb('R','T','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_sgemlqt
     !> DGEMLQT: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'T':   Q**T C            C Q**T
     !> where Q is a real orthogonal matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**T
     !> generated using the compact WY representation as returned by DGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_dgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           real(dp),intent(in) :: v(ldv,*),t(ldt,*)
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_dlarfb('L','T','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_dlarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_dlarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_dlarfb('R','T','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_dgemlqt
#ifdef LA_WITH_XDP
     !> XGEMLQT: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'T':   Q**T C            C Q**T
     !> where Q is a real orthogonal matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**T
     !> generated using the compact WY representation as returned by XGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_xgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           real(xdp),intent(in) :: v(ldv,*),t(ldt,*)
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('XGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_xlarfb('L','T','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_xlarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_xlarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_xlarfb('R','T','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_xgemlqt
#endif
#ifdef LA_WITH_QP
     !> QGEMLQT: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'T':   Q**T C            C Q**T
     !> where Q is a real orthogonal matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**T
     !> generated using the compact WY representation as returned by QGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_qgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           real(qp),intent(in) :: v(ldv,*),t(ldt,*)
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_qlarfb('L','T','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_qlarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_qlarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_qlarfb('R','T','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_qgemlqt
#endif

     !> STPLQT2: computes a LQ a factorization of a real "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_stplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           real(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('STPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_slarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              if (i < m) then
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_sgemv('N',m - i,p,one,b(i + 1,1),ldb,b(i,1),ldb,one,t(m, &
                           1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_sger(m - i,p,alpha,t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n) * (alpha * c(i,i:n)^h)
              alpha = -t(1,i)
              do j = 1,i - 1
                 t(i,j) = zero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = alpha*b(i,n - l + j)
              end do
              call la_strmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_sgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,zero,t( &
                         i,mp),ldt)
              ! b1
              call la_sgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,one,t(i,1),ldt &
                        )
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
             call la_strmv('L','T','N',i - 1,t,ldt,t(i,1),ldt)
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = zero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = t(j,i)
                 t(j,i) = zero
              end do
           end do
     end subroutine la_stplqt2
     !> DTPLQT2: computes a LQ a factorization of a real "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_dtplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           real(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_dlarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              if (i < m) then
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_dgemv('N',m - i,p,one,b(i + 1,1),ldb,b(i,1),ldb,one,t(m, &
                           1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_dger(m - i,p,alpha,t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n) * (alpha * c(i,i:n)^h)
              alpha = -t(1,i)
              do j = 1,i - 1
                 t(i,j) = zero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = alpha*b(i,n - l + j)
              end do
              call la_dtrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_dgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,zero,t( &
                         i,mp),ldt)
              ! b1
              call la_dgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,one,t(i,1),ldt &
                        )
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
             call la_dtrmv('L','T','N',i - 1,t,ldt,t(i,1),ldt)
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = zero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = t(j,i)
                 t(j,i) = zero
              end do
           end do
     end subroutine la_dtplqt2
#ifdef LA_WITH_XDP
     !> XTPLQT2: computes a LQ a factorization of a real "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_xtplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           real(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_xlarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              if (i < m) then
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_xgemv('N',m - i,p,one,b(i + 1,1),ldb,b(i,1),ldb,one,t(m, &
                           1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_xger(m - i,p,alpha,t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n) * (alpha * c(i,i:n)^h)
              alpha = -t(1,i)
              do j = 1,i - 1
                 t(i,j) = zero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = alpha*b(i,n - l + j)
              end do
              call la_xtrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_xgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,zero,t( &
                         i,mp),ldt)
              ! b1
              call la_xgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,one,t(i,1),ldt &
                        )
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
             call la_xtrmv('L','T','N',i - 1,t,ldt,t(i,1),ldt)
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = zero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = t(j,i)
                 t(j,i) = zero
              end do
           end do
     end subroutine la_xtplqt2
#endif
#ifdef LA_WITH_QP
     !> QTPLQT2: computes a LQ a factorization of a real "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_qtplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           real(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_qlarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              if (i < m) then
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_qgemv('N',m - i,p,one,b(i + 1,1),ldb,b(i,1),ldb,one,t(m, &
                           1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_qger(m - i,p,alpha,t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n) * (alpha * c(i,i:n)^h)
              alpha = -t(1,i)
              do j = 1,i - 1
                 t(i,j) = zero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = alpha*b(i,n - l + j)
              end do
              call la_qtrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_qgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,zero,t( &
                         i,mp),ldt)
              ! b1
              call la_qgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,one,t(i,1),ldt &
                        )
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
             call la_qtrmv('L','T','N',i - 1,t,ldt,t(i,1),ldt)
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = zero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = t(j,i)
                 t(j,i) = zero
              end do
           end do
     end subroutine la_qtplqt2
#endif

     !> STPMLQT: applies a real orthogonal matrix Q obtained from a
     !> "triangular-pentagonal" real block reflector H to a general
     !> real matrix C, which consists of two blocks A and B.

     pure subroutine la_stpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           real(sp),intent(in) :: v(ldv,*),t(ldt,*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('STPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_stprfb('L','T','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_stprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_stprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_stprfb('R','T','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_stpmlqt
     !> DTPMQRT applies a real orthogonal matrix Q obtained from a
     !> "triangular-pentagonal" real block reflector H to a general
     !> real matrix C, which consists of two blocks A and B.

     pure subroutine la_dtpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           real(dp),intent(in) :: v(ldv,*),t(ldt,*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('DTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_dtprfb('L','T','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_dtprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_dtprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_dtprfb('R','T','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_dtpmlqt
#ifdef LA_WITH_XDP
     !> XTPMQRT applies a real orthogonal matrix Q obtained from a
     !> "triangular-pentagonal" real block reflector H to a general
     !> real matrix C, which consists of two blocks A and B.

     pure subroutine la_xtpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           real(xdp),intent(in) :: v(ldv,*),t(ldt,*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('XTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_xtprfb('L','T','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_xtprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_xtprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_xtprfb('R','T','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_xtpmlqt
#endif
#ifdef LA_WITH_QP
     !> QTPMQRT applies a real orthogonal matrix Q obtained from a
     !> "triangular-pentagonal" real block reflector H to a general
     !> real matrix C, which consists of two blocks A and B.

     pure subroutine la_qtpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           real(qp),intent(in) :: v(ldv,*),t(ldt,*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'T')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('QTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_qtprfb('L','T','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_qtprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_qtprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_qtprfb('R','T','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_qtpmlqt
#endif

     !> SGELQ2: computes an LQ factorization of a real m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_sgelq2(m,n,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_slarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 aii = a(i,i)
                 a(i,i) = one
                 call la_slarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
                 a(i,i) = aii
              end if
           end do
           return
     end subroutine la_sgelq2
     !> DGELQ2: computes an LQ factorization of a real m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_dgelq2(m,n,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_dlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 aii = a(i,i)
                 a(i,i) = one
                 call la_dlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
                 a(i,i) = aii
              end if
           end do
           return
     end subroutine la_dgelq2
#ifdef LA_WITH_XDP
     !> XGELQ2: computes an LQ factorization of a real m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_xgelq2(m,n,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(xdp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('XGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_xlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 aii = a(i,i)
                 a(i,i) = one
                 call la_xlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
                 a(i,i) = aii
              end if
           end do
           return
     end subroutine la_xgelq2
#endif
#ifdef LA_WITH_QP
     !> QGELQ2: computes an LQ factorization of a real m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_qgelq2(m,n,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_qlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 aii = a(i,i)
                 a(i,i) = one
                 call la_qlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
                 a(i,i) = aii
              end if
           end do
           return
     end subroutine la_qgelq2
#endif

     !> SGELQF: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_sgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'SGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'SGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'SGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_sgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_slarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_slarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_sgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_sgelqf
     !> DGELQF: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_dgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'DGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'DGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'DGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_dgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_dlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_dlarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_dgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_dgelqf
#ifdef LA_WITH_XDP
     !> XGELQF: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_xgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'XGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'XGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'XGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_xgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_xlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_xlarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_xgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_xgelqf
#endif
#ifdef LA_WITH_QP
     !> QGELQF: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_qgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'QGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'QGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'QGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_qgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_qlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_qlarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_qgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_qgelqf
#endif

     !> SGELQT3: recursively computes a LQ factorization of a real M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_sgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_slarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_sgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) = q1^h a(j1:m,1:n) [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_strmm('R','U','T','U',m2,m1,one,a,lda,t(i1,1),ldt)
              call la_sgemm('N','T',m2,m1,n - m1,one,a(i1,i1),lda,a(1,i1),lda, &
                        one,t(i1,1),ldt)
              call la_strmm('R','U','N','N',m2,m1,one,t,ldt,t(i1,1),ldt)
              call la_sgemm('N','N',m2,n - m1,m1,-one,t(i1,1),ldt,a(1,i1),lda, &
                        one,a(i1,i1),lda)
              call la_strmm('R','U','N','U',m2,m1,one,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = 0
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_sgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_strmm('R','U','T','U',m1,m2,one,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_sgemm('N','T',m1,m2,n - m,one,a(1,j1),lda,a(i1,j1),lda, &
                        one,t(1,i1),ldt)
              call la_strmm('L','U','N','N',m1,m2,-one,t,ldt,t(1,i1),ldt)

              call la_strmm('R','U','N','N',m1,m2,one,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_sgelqt3
     !> DGELQT3: recursively computes a LQ factorization of a real M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_dgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_dlarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_dgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) = q1^h a(j1:m,1:n) [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_dtrmm('R','U','T','U',m2,m1,one,a,lda,t(i1,1),ldt)
              call la_dgemm('N','T',m2,m1,n - m1,one,a(i1,i1),lda,a(1,i1),lda, &
                        one,t(i1,1),ldt)
              call la_dtrmm('R','U','N','N',m2,m1,one,t,ldt,t(i1,1),ldt)
              call la_dgemm('N','N',m2,n - m1,m1,-one,t(i1,1),ldt,a(1,i1),lda, &
                        one,a(i1,i1),lda)
              call la_dtrmm('R','U','N','U',m2,m1,one,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = 0
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_dgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_dtrmm('R','U','T','U',m1,m2,one,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_dgemm('N','T',m1,m2,n - m,one,a(1,j1),lda,a(i1,j1),lda, &
                        one,t(1,i1),ldt)
              call la_dtrmm('L','U','N','N',m1,m2,-one,t,ldt,t(1,i1),ldt)

              call la_dtrmm('R','U','N','N',m1,m2,one,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_dgelqt3
#ifdef LA_WITH_XDP
     !> XGELQT3: recursively computes a LQ factorization of a real M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_xgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('XGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_xlarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_xgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) = q1^h a(j1:m,1:n) [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_xtrmm('R','U','T','U',m2,m1,one,a,lda,t(i1,1),ldt)
              call la_xgemm('N','T',m2,m1,n - m1,one,a(i1,i1),lda,a(1,i1),lda, &
                        one,t(i1,1),ldt)
              call la_xtrmm('R','U','N','N',m2,m1,one,t,ldt,t(i1,1),ldt)
              call la_xgemm('N','N',m2,n - m1,m1,-one,t(i1,1),ldt,a(1,i1),lda, &
                        one,a(i1,i1),lda)
              call la_xtrmm('R','U','N','U',m2,m1,one,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = 0
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_xgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_xtrmm('R','U','T','U',m1,m2,one,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_xgemm('N','T',m1,m2,n - m,one,a(1,j1),lda,a(i1,j1),lda, &
                        one,t(1,i1),ldt)
              call la_xtrmm('L','U','N','N',m1,m2,-one,t,ldt,t(1,i1),ldt)

              call la_xtrmm('R','U','N','N',m1,m2,one,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_xgelqt3
#endif
#ifdef LA_WITH_QP
     !> QGELQT3: recursively computes a LQ factorization of a real M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_qgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_qlarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_qgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) = q1^h a(j1:m,1:n) [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_qtrmm('R','U','T','U',m2,m1,one,a,lda,t(i1,1),ldt)
              call la_qgemm('N','T',m2,m1,n - m1,one,a(i1,i1),lda,a(1,i1),lda, &
                        one,t(i1,1),ldt)
              call la_qtrmm('R','U','N','N',m2,m1,one,t,ldt,t(i1,1),ldt)
              call la_qgemm('N','N',m2,n - m1,m1,-one,t(i1,1),ldt,a(1,i1),lda, &
                        one,a(i1,i1),lda)
              call la_qtrmm('R','U','N','U',m2,m1,one,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = 0
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_qgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_qtrmm('R','U','T','U',m1,m2,one,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_qgemm('N','T',m1,m2,n - m,one,a(1,j1),lda,a(i1,j1),lda, &
                        one,t(1,i1),ldt)
              call la_qtrmm('L','U','N','N',m1,m2,-one,t,ldt,t(1,i1),ldt)

              call la_qtrmm('R','U','N','N',m1,m2,one,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_qgelqt3
#endif

     !> SGEQL2: computes a QL factorization of a real m by n matrix A:
     !> A = Q * L.

     pure subroutine la_sgeql2(m,n,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              call la_slarfg(m - k + i,a(m - k + i,n - k + i),a(1,n - k + i),1,tau(i))
              ! apply h(i) to a(1:m-k+i,1:n-k+i-1) from the left
              aii = a(m - k + i,n - k + i)
              a(m - k + i,n - k + i) = one
              call la_slarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,tau(i),a,lda,work)

              a(m - k + i,n - k + i) = aii
           end do
           return
     end subroutine la_sgeql2
     !> DGEQL2: computes a QL factorization of a real m by n matrix A:
     !> A = Q * L.

     pure subroutine la_dgeql2(m,n,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              call la_dlarfg(m - k + i,a(m - k + i,n - k + i),a(1,n - k + i),1,tau(i))
              ! apply h(i) to a(1:m-k+i,1:n-k+i-1) from the left
              aii = a(m - k + i,n - k + i)
              a(m - k + i,n - k + i) = one
              call la_dlarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,tau(i),a,lda,work)

              a(m - k + i,n - k + i) = aii
           end do
           return
     end subroutine la_dgeql2
#ifdef LA_WITH_XDP
     !> XGEQL2: computes a QL factorization of a real m by n matrix A:
     !> A = Q * L.

     pure subroutine la_xgeql2(m,n,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(xdp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('XGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              call la_xlarfg(m - k + i,a(m - k + i,n - k + i),a(1,n - k + i),1,tau(i))
              ! apply h(i) to a(1:m-k+i,1:n-k+i-1) from the left
              aii = a(m - k + i,n - k + i)
              a(m - k + i,n - k + i) = one
              call la_xlarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,tau(i),a,lda,work)

              a(m - k + i,n - k + i) = aii
           end do
           return
     end subroutine la_xgeql2
#endif
#ifdef LA_WITH_QP
     !> QGEQL2: computes a QL factorization of a real m by n matrix A:
     !> A = Q * L.

     pure subroutine la_qgeql2(m,n,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              call la_qlarfg(m - k + i,a(m - k + i,n - k + i),a(1,n - k + i),1,tau(i))
              ! apply h(i) to a(1:m-k+i,1:n-k+i-1) from the left
              aii = a(m - k + i,n - k + i)
              a(m - k + i,n - k + i) = one
              call la_qlarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,tau(i),a,lda,work)

              a(m - k + i,n - k + i) = aii
           end do
           return
     end subroutine la_qgeql2
#endif

     !> SGEQLF: computes a QL factorization of a real M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_sgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'SGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'SGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'SGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_sgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_slarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**t to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_slarfb('LEFT','TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - 1, &
                    n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_sgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_sgeqlf
     !> DGEQLF: computes a QL factorization of a real M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_dgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'DGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'DGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'DGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_dgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_dlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**t to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_dlarfb('LEFT','TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - 1, &
                    n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_dgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_dgeqlf
#ifdef LA_WITH_XDP
     !> XGEQLF: computes a QL factorization of a real M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_xgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'XGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'XGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'XGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_xgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_xlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**t to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_xlarfb('LEFT','TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - 1, &
                    n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_xgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_xgeqlf
#endif
#ifdef LA_WITH_QP
     !> QGEQLF: computes a QL factorization of a real M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_qgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'QGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'QGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'QGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_qgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_qlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**t to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_qlarfb('LEFT','TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - 1, &
                    n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_qgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_qgeqlf
#endif

     !> SLAMSWLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (SLASWLQ)

     pure subroutine la_slamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),t(ldt,*)
           real(sp),intent(out) :: work(*)
           real(sp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,lw,ctr
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('SLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_sgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_stpmlqt('L','T',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_stpmlqt('L','T',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                           1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_sgemlqt('L','T',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_sgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_stpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c( &
                         1,1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_stpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_stpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
                  ctr = ctr - 1
                  call la_stpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                            c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_sgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ii = n - kk + 1
              ctr = 1
              call la_sgemlqt('R','T',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_stpmlqt('R','T',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_stpmlqt('R','T',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1,1), &
                          ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_slamswlq
     !> DLAMSWLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (DLASWLQ)

     pure subroutine la_dlamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),t(ldt,*)
           real(dp),intent(out) :: work(*)
           real(dp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,ctr,lw
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('DLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_dgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_dtpmlqt('L','T',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_dtpmlqt('L','T',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_dgemlqt('L','T',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_dgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_dtpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_dtpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_dtpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                           c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
                  ctr = ctr - 1
                  call la_dtpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                            c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_dgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ctr = 1
              ii = n - kk + 1
              call la_dgemlqt('R','T',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_dtpmlqt('R','T',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_dtpmlqt('R','T',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1,1), &
                          ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_dlamswlq
#ifdef LA_WITH_XDP
     !> XLAMSWLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (XLASWLQ)

     pure subroutine la_xlamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),t(ldt,*)
           real(xdp),intent(out) :: work(*)
           real(xdp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,ctr,lw
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('XLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_xgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_xtpmlqt('L','T',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_xtpmlqt('L','T',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_xgemlqt('L','T',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_xgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_xtpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_xtpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_xtpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                           c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
                  ctr = ctr - 1
                  call la_xtpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                            c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_xgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ctr = 1
              ii = n - kk + 1
              call la_xgemlqt('R','T',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_xtpmlqt('R','T',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_xtpmlqt('R','T',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1,1), &
                          ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_xlamswlq
#endif
#ifdef LA_WITH_QP
     !> QLAMSWLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (QLASWLQ)

     pure subroutine la_qlamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),t(ldt,*)
           real(qp),intent(out) :: work(*)
           real(qp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,ctr,lw
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('QLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_qgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_qtpmlqt('L','T',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_qtpmlqt('L','T',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_qgemlqt('L','T',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_qgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_qtpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_qtpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_qtpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                           c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
                  ctr = ctr - 1
                  call la_qtpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                            c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_qgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ctr = 1
              ii = n - kk + 1
              call la_qgemlqt('R','T',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_qtpmlqt('R','T',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_qtpmlqt('R','T',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1,1), &
                          ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_qlamswlq
#endif

     !> STPLQT: computes a blocked LQ factorization of a real
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_stplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('STPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_stplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_stprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_stplqt
     !> DTPLQT: computes a blocked LQ factorization of a real
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_dtplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('DTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_dtplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_dtprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_dtplqt
#ifdef LA_WITH_XDP
     !> XTPLQT: computes a blocked LQ factorization of a real
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_xtplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('XTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_xtplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_xtprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_xtplqt
#endif
#ifdef LA_WITH_QP
     !> QTPLQT: computes a blocked LQ factorization of a real
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_qtplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('QTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_qtplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_qtprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_qtplqt
#endif

     !> DGELQT computes a blocked LQ factorization of a real M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_sgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_sgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_slarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_sgelqt
     !> DGELQT: computes a blocked LQ factorization of a real M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_dgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_dgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_dlarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_dgelqt
#ifdef LA_WITH_XDP
     !> XGELQT: computes a blocked LQ factorization of a real M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_xgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_xgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_xlarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_xgelqt
#endif
#ifdef LA_WITH_QP
     !> QGELQT: computes a blocked LQ factorization of a real M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_qgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_qgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_qlarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_qgelqt
#endif

     !> SGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product
     !> of blocked elementary reflectors computed by short wide LQ
     !> factorization (SGELQ)

     pure subroutine la_sgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),t(*)
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = real(lw,KIND=sp)
           end if
           if (info /= 0) then
             call la_xerbla('SGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_sgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_slamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = real(lw,KIND=sp)
           return
     end subroutine la_sgemlq
     !> DGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product
     !> of blocked elementary reflectors computed by short wide LQ
     !> factorization (DGELQ)

     pure subroutine la_dgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),t(*)
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('DGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_dgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_dlamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_dgemlq
#ifdef LA_WITH_XDP
     !> XGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product
     !> of blocked elementary reflectors computed by short wide LQ
     !> factorization (XGELQ)

     pure subroutine la_xgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),t(*)
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('XGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_xgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_xlamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_xgemlq
#endif
#ifdef LA_WITH_QP
     !> QGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product
     !> of blocked elementary reflectors computed by short wide LQ
     !> factorization (QGELQ)

     pure subroutine la_qgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),t(*)
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'T')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('QGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_qgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_qlamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_qgemlq
#endif

     !> SLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a real M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_slaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb <= 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('SLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_sgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_sgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_stplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_stplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_slaswlq
     !> DLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a real M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_dlaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb < 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('DLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_dgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_dgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_dtplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_dtplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_dlaswlq
#ifdef LA_WITH_XDP
     !> XLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a real M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_xlaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb < 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('XLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_xgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_xgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_xtplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_xtplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_xlaswlq
#endif
#ifdef LA_WITH_QP
     !> QLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a real M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_qlaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb < 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('QLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_qgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_qgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_qtplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_qtplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_qlaswlq
#endif

     !> SGELQ: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_sgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'SGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'SGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('SGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_sgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_slaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_sgelq
     !> DGELQ: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_dgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'DGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'DGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('DGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_dgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_dlaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_dgelq
#ifdef LA_WITH_XDP
     !> XGELQ: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_xgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'XGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'XGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('XGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_xgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_xlaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_xgelq
#endif
#ifdef LA_WITH_QP
     !> QGELQ: computes an LQ factorization of a real M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_qgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'QGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'QGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('QGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_qgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_qlaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_qgelq
#endif

     !> CTPLQT2: computes a LQ a factorization of a complex "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_ctplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_clarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              t(1,i) = conjg(t(1,i))
              if (i < m) then
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_cgemv('N',m - i,p,cone,b(i + 1,1),ldb,b(i,1),ldb,cone,t( &
                           m,1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_cgerc(m - i,p, (alpha),t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n)**h * (alpha * c(i,i:n))
              alpha = -(t(1,i))
              do j = 1,i - 1
                 t(i,j) = czero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              do j = 1,n - l + p
                b(i,j) = conjg(b(i,j))
              end do
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = (alpha*b(i,n - l + j))
              end do
              call la_ctrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_cgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,czero, &
                        t(i,mp),ldt)
              ! b1
              call la_cgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,cone,t(i,1), &
                        ldt)
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              call la_ctrmv('L','C','N',i - 1,t,ldt,t(i,1),ldt)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              do j = 1,n - l + p
                 b(i,j) = conjg(b(i,j))
              end do
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = czero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = (t(j,i))
                 t(j,i) = czero
              end do
           end do
     end subroutine la_ctplqt2
     !> ZTPLQT2: computes a LQ a factorization of a complex "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_ztplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_zlarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              t(1,i) = conjg(t(1,i))
              if (i < m) then
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_zgemv('N',m - i,p,cone,b(i + 1,1),ldb,b(i,1),ldb,cone,t( &
                           m,1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_zgerc(m - i,p, (alpha),t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n)**h * (alpha * c(i,i:n))
              alpha = -(t(1,i))
              do j = 1,i - 1
                 t(i,j) = czero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              do j = 1,n - l + p
                b(i,j) = conjg(b(i,j))
              end do
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = (alpha*b(i,n - l + j))
              end do
              call la_ztrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_zgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,czero, &
                        t(i,mp),ldt)
              ! b1
              call la_zgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,cone,t(i,1), &
                        ldt)
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              call la_ztrmv('L','C','N',i - 1,t,ldt,t(i,1),ldt)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              do j = 1,n - l + p
                 b(i,j) = conjg(b(i,j))
              end do
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = czero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = (t(j,i))
                 t(j,i) = czero
              end do
           end do
     end subroutine la_ztplqt2
#ifdef LA_WITH_XDP
     !> YTPLQT2: computes a LQ a factorization of a complex "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_ytplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           complex(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_ylarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              t(1,i) = conjg(t(1,i))
              if (i < m) then
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_ygemv('N',m - i,p,cone,b(i + 1,1),ldb,b(i,1),ldb,cone,t( &
                           m,1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_ygerc(m - i,p, (alpha),t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n)**h * (alpha * c(i,i:n))
              alpha = -(t(1,i))
              do j = 1,i - 1
                 t(i,j) = czero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              do j = 1,n - l + p
                b(i,j) = conjg(b(i,j))
              end do
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = (alpha*b(i,n - l + j))
              end do
              call la_ytrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_ygemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,czero, &
                        t(i,mp),ldt)
              ! b1
              call la_ygemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,cone,t(i,1), &
                        ldt)
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              call la_ytrmv('L','C','N',i - 1,t,ldt,t(i,1),ldt)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              do j = 1,n - l + p
                 b(i,j) = conjg(b(i,j))
              end do
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = czero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = (t(j,i))
                 t(j,i) = czero
              end do
           end do
     end subroutine la_ytplqt2
#endif
#ifdef LA_WITH_QP
     !> WTPLQT2: computes a LQ a factorization of a complex "triangular-pentagonal"
     !> matrix C, which is composed of a triangular block A and pentagonal block B,
     !> using the compact WY representation for Q.

     pure subroutine la_wtplqt2(m,n,l,a,lda,b,ldb,t,ldt,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,p,mp,np
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. l > min(m,n)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldb < max(1,m)) then
              info = -7
           else if (ldt < max(1,m)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WTPLQT2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           do i = 1,m
              ! generate elementary reflector h(i) to annihilate b(i,:)
              p = n - l + min(l,i)
              call la_wlarfg(p + 1,a(i,i),b(i,1),ldb,t(1,i))
              t(1,i) = conjg(t(1,i))
              if (i < m) then
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
                 ! w(m-i:1) := c(i+1:m,i:n) * c(i,i:n) [use w = t(m,:)]
                 do j = 1,m - i
                    t(m,j) = (a(i + j,i))
                 end do
                 call la_wgemv('N',m - i,p,cone,b(i + 1,1),ldb,b(i,1),ldb,cone,t( &
                           m,1),ldt)
                 ! c(i+1:m,i:n) = c(i+1:m,i:n) + alpha * c(i,i:n)*w(m-1:1)^h
                 alpha = -(t(1,i))
                 do j = 1,m - i
                    a(i + j,i) = a(i + j,i) + alpha*(t(m,j))
                 end do
                 call la_wgerc(m - i,p, (alpha),t(m,1),ldt,b(i,1),ldb,b(i + 1,1), &
                           ldb)
                 do j = 1,p
                    b(i,j) = conjg(b(i,j))
                 end do
              end if
           end do
           do i = 2,m
              ! t(i,1:i-1) := c(i:i-1,1:n)**h * (alpha * c(i,i:n))
              alpha = -(t(1,i))
              do j = 1,i - 1
                 t(i,j) = czero
              end do
              p = min(i - 1,l)
              np = min(n - l + 1,n)
              mp = min(p + 1,m)
              do j = 1,n - l + p
                b(i,j) = conjg(b(i,j))
              end do
              ! triangular part of b2
              do j = 1,p
                 t(i,j) = (alpha*b(i,n - l + j))
              end do
              call la_wtrmv('L','N','N',p,b(1,np),ldb,t(i,1),ldt)
              ! rectangular part of b2
              call la_wgemv('N',i - 1 - p,l,alpha,b(mp,np),ldb,b(i,np),ldb,czero, &
                        t(i,mp),ldt)
              ! b1
              call la_wgemv('N',i - 1,n - l,alpha,b,ldb,b(i,1),ldb,cone,t(i,1), &
                        ldt)
              ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(i,1:i-1)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              call la_wtrmv('L','C','N',i - 1,t,ldt,t(i,1),ldt)
              do j = 1,i - 1
                 t(i,j) = conjg(t(i,j))
              end do
              do j = 1,n - l + p
                 b(i,j) = conjg(b(i,j))
              end do
              ! t(i,i) = tau(i)
              t(i,i) = t(1,i)
              t(1,i) = czero
           end do
           do i = 1,m
              do j = i + 1,m
                 t(i,j) = (t(j,i))
                 t(j,i) = czero
              end do
           end do
     end subroutine la_wtplqt2
#endif

     !> CUNG2L: generates an m by n complex matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by CGEQLF.

     pure subroutine la_cung2l(m,n,k,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CUNG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = czero
              end do
              a(m - n + j,j) = cone
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = cone
              call la_clarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_cscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = cone - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to czero
              do l = m - n + ii + 1,m
                 a(l,ii) = czero
              end do
           end do
           return
     end subroutine la_cung2l
     !> ZUNG2L: generates an m by n complex matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by ZGEQLF.

     pure subroutine la_zung2l(m,n,k,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZUNG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = czero
              end do
              a(m - n + j,j) = cone
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = cone
              call la_zlarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_zscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = cone - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to czero
              do l = m - n + ii + 1,m
                 a(l,ii) = czero
              end do
           end do
           return
     end subroutine la_zung2l
#ifdef LA_WITH_XDP
     !> YUNG2L: generates an m by n complex matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by YGEQLF.

     pure subroutine la_yung2l(m,n,k,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('YUNG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = czero
              end do
              a(m - n + j,j) = cone
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = cone
              call la_ylarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_yscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = cone - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to czero
              do l = m - n + ii + 1,m
                 a(l,ii) = czero
              end do
           end do
           return
     end subroutine la_yung2l
#endif
#ifdef LA_WITH_QP
     !> WUNG2L: generates an m by n complex matrix Q with orthonormal columns,
     !> which is defined as the last n columns of a product of k elementary
     !> reflectors of order m
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by WGEQLF.

     pure subroutine la_wung2l(m,n,k,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,j,l
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WUNG2L',-info)
              return
           end if
           ! quick return if possible
           if (n <= 0) return
           ! initialise columns 1:n-k to columns of the unit matrix
           do j = 1,n - k
              do l = 1,m
                 a(l,j) = czero
              end do
              a(m - n + j,j) = cone
           end do
           do i = 1,k
              ii = n - k + i
              ! apply h(i) to a(1:m-k+i,1:n-k+i) from the left
              a(m - n + ii,ii) = cone
              call la_wlarf('LEFT',m - n + ii,ii - 1,a(1,ii),1,tau(i),a,lda,work)

              call la_wscal(m - n + ii - 1,-tau(i),a(1,ii),1)
              a(m - n + ii,ii) = cone - tau(i)
              ! set a(m-k+i+1:m,n-k+i) to czero
              do l = m - n + ii + 1,m
                 a(l,ii) = czero
              end do
           end do
           return
     end subroutine la_wung2l
#endif

     !> CUNGL2: generates an m-by-n complex matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by CGELQF.

     pure subroutine la_cungl2(m,n,k,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CUNGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = czero
                 end do
                 if (j > k .and. j <= m) a(j,j) = cone
              end do
           end if
           do i = k,1,-1
              ! apply h(i)**h to a(i:m,i:n) from the right
              if (i < n) then
                 call la_clacgv(n - i,a(i,i + 1),lda)
                 if (i < m) then
                    a(i,i) = cone
                    call la_clarf('RIGHT',m - i,n - i + 1,a(i,i),lda,conjg(tau(i)),a( &
                              i + 1,i),lda,work)
                 end if
                 call la_cscal(n - i,-tau(i),a(i,i + 1),lda)
                 call la_clacgv(n - i,a(i,i + 1),lda)
              end if
              a(i,i) = cone - conjg(tau(i))
              ! set a(i,1:i-1,i) to czero
              do l = 1,i - 1
                 a(i,l) = czero
              end do
           end do
           return
     end subroutine la_cungl2
     !> ZUNGL2: generates an m-by-n complex matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by ZGELQF.

     pure subroutine la_zungl2(m,n,k,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZUNGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = czero
                 end do
                 if (j > k .and. j <= m) a(j,j) = cone
              end do
           end if
           do i = k,1,-1
              ! apply h(i)**h to a(i:m,i:n) from the right
              if (i < n) then
                 call la_zlacgv(n - i,a(i,i + 1),lda)
                 if (i < m) then
                    a(i,i) = cone
                    call la_zlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,conjg(tau(i)),a( &
                              i + 1,i),lda,work)
                 end if
                 call la_zscal(n - i,-tau(i),a(i,i + 1),lda)
                 call la_zlacgv(n - i,a(i,i + 1),lda)
              end if
              a(i,i) = cone - conjg(tau(i))
              ! set a(i,1:i-1) to czero
              do l = 1,i - 1
                 a(i,l) = czero
              end do
           end do
           return
     end subroutine la_zungl2
#ifdef LA_WITH_XDP
     !> YUNGL2: generates an m-by-n complex matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by YGELQF.

     pure subroutine la_yungl2(m,n,k,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('YUNGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = czero
                 end do
                 if (j > k .and. j <= m) a(j,j) = cone
              end do
           end if
           do i = k,1,-1
              ! apply h(i)**h to a(i:m,i:n) from the right
              if (i < n) then
                 call la_ylacgv(n - i,a(i,i + 1),lda)
                 if (i < m) then
                    a(i,i) = cone
                    call la_ylarf('RIGHT',m - i,n - i + 1,a(i,i),lda,conjg(tau(i)),a( &
                              i + 1,i),lda,work)
                 end if
                 call la_yscal(n - i,-tau(i),a(i,i + 1),lda)
                 call la_ylacgv(n - i,a(i,i + 1),lda)
              end if
              a(i,i) = cone - conjg(tau(i))
              ! set a(i,1:i-1) to czero
              do l = 1,i - 1
                 a(i,l) = czero
              end do
           end do
           return
     end subroutine la_yungl2
#endif
#ifdef LA_WITH_QP
     !> WUNGL2: generates an m-by-n complex matrix Q with orthonormal rows,
     !> which is defined as the first m rows of a product of k elementary
     !> reflectors of order n
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by WGELQF.

     pure subroutine la_wungl2(m,n,k,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WUNGL2',-info)
              return
           end if
           ! quick return if possible
           if (m <= 0) return
           if (k < m) then
              ! initialise rows k+1:m to rows of the unit matrix
              do j = 1,n
                 do l = k + 1,m
                    a(l,j) = czero
                 end do
                 if (j > k .and. j <= m) a(j,j) = cone
              end do
           end if
           do i = k,1,-1
              ! apply h(i)**h to a(i:m,i:n) from the right
              if (i < n) then
                 call la_wlacgv(n - i,a(i,i + 1),lda)
                 if (i < m) then
                    a(i,i) = cone
                    call la_wlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,conjg(tau(i)),a( &
                              i + 1,i),lda,work)
                 end if
                 call la_wscal(n - i,-tau(i),a(i,i + 1),lda)
                 call la_wlacgv(n - i,a(i,i + 1),lda)
              end if
              a(i,i) = cone - conjg(tau(i))
              ! set a(i,1:i-1) to czero
              do l = 1,i - 1
                 a(i,l) = czero
              end do
           end do
           return
     end subroutine la_wungl2
#endif

     !> CUNGLQ: generates an M-by-N complex matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by CGELQF.

     pure subroutine la_cunglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'CUNGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CUNGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'CUNGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'CUNGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to czero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_cungl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_clarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**h to a(i+ib:m,i:n) from the right
                    call la_clarfb('RIGHT','CONJUGATE TRANSPOSE','FORWARD','ROWWISE',m - i - &
                    ib + 1,n - i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1), &
                              ldwork)
                 end if
                 ! apply h**h to columns i:n of current block
                 call la_cungl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to czero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_cunglq
     !> ZUNGLQ: generates an M-by-N complex matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by ZGELQF.

     pure subroutine la_zunglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'ZUNGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZUNGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'ZUNGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'ZUNGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to czero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_zungl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_zlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**h to a(i+ib:m,i:n) from the right
                    call la_zlarfb('RIGHT','CONJUGATE TRANSPOSE','FORWARD','ROWWISE',m - i - &
                    ib + 1,n - i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1), &
                              ldwork)
                 end if
                 ! apply h**h to columns i:n of current block
                 call la_zungl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to czero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_zunglq
#ifdef LA_WITH_XDP
     !> YUNGLQ: generates an M-by-N complex matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by YGELQF.

     pure subroutine la_yunglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'YUNGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('YUNGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'YUNGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'YUNGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to czero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_yungl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_ylarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**h to a(i+ib:m,i:n) from the right
                    call la_ylarfb('RIGHT','CONJUGATE TRANSPOSE','FORWARD','ROWWISE',m - i - &
                    ib + 1,n - i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1), &
                              ldwork)
                 end if
                 ! apply h**h to columns i:n of current block
                 call la_yungl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to czero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_yunglq
#endif
#ifdef LA_WITH_QP
     !> WUNGLQ: generates an M-by-N complex matrix Q with orthonormal rows,
     !> which is defined as the first M rows of a product of K elementary
     !> reflectors of order N
     !> Q  =  H(k)**H . . . H(2)**H H(1)**H
     !> as returned by WGELQF.

     pure subroutine la_wunglq(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,ki,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'WUNGLQ',' ',m,n,k,-1)
           lwkopt = max(1,m)*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (k < 0 .or. k > m) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WUNGLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m <= 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'WUNGLQ',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'WUNGLQ',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the last block.
              ! the first kk rows are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              ! set a(kk+1:m,1:kk) to czero.
              do j = 1,kk
                 do i = kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the last or only block.
           if (kk < m) call la_wungl2(m - kk,n - kk,k - kk,a(kk + 1,kk + 1),lda,tau(kk + 1),work, &
                      iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = ki + 1,1,-nb
                 ib = min(nb,k - i + 1)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_wlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h**h to a(i+ib:m,i:n) from the right
                    call la_wlarfb('RIGHT','CONJUGATE TRANSPOSE','FORWARD','ROWWISE',m - i - &
                    ib + 1,n - i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1), &
                              ldwork)
                 end if
                 ! apply h**h to columns i:n of current block
                 call la_wungl2(ib,n - i + 1,ib,a(i,i),lda,tau(i),work,iinfo)
                 ! set columns 1:i-1 of current block to czero
                 do j = 1,i - 1
                    do l = i,i + ib - 1
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_wunglq
#endif

     !> CUNGQL: generates an M-by-N complex matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by CGEQLF.

     pure subroutine la_cungql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'CUNGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CUNGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'CUNGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'CUNGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to czero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_cung2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_clarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_clarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_cung2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to czero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_cungql
     !> ZUNGQL: generates an M-by-N complex matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by ZGEQLF.

     pure subroutine la_zungql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'ZUNGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZUNGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'ZUNGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'ZUNGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to czero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_zung2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_zlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_zlarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_zung2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to czero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_zungql
#ifdef LA_WITH_XDP
     !> YUNGQL: generates an M-by-N complex matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by YGEQLF.

     pure subroutine la_yungql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'YUNGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YUNGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'YUNGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'YUNGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to czero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_yung2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_ylarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_ylarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_yung2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to czero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_yungql
#endif
#ifdef LA_WITH_QP
     !> WUNGQL: generates an M-by-N complex matrix Q with orthonormal columns,
     !> which is defined as the last N columns of a product of K elementary
     !> reflectors of order M
     !> Q  =  H(k) . . . H(2) H(1)
     !> as returned by WGEQLF.

     pure subroutine la_wungql(m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,j,kk,l,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0 .or. n > m) then
              info = -2
           else if (k < 0 .or. k > n) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'WUNGQL',' ',m,n,k,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -8
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WUNGQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           nbmin = 2
           nx = 0
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'WUNGQL',' ',m,n,k,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'WUNGQL',' ',m,n,k,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code after the first block.
              ! the last kk columns are handled by the block method.
              kk = min(k, ((k - nx + nb - 1)/nb)*nb)
              ! set a(m-kk+1:m,1:n-kk) to czero.
              do j = 1,n - kk
                 do i = m - kk + 1,m
                    a(i,j) = czero
                 end do
              end do
           else
              kk = 0
           end if
           ! use unblocked code for the first or only block.
           call la_wung2l(m - kk,n - kk,k - kk,a,lda,tau,work,iinfo)
           if (kk > 0) then
              ! use blocked code
              do i = k - kk + 1,k,nb
                 ib = min(nb,k - i + 1)
                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_wlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_wlarfb('LEFT','NO TRANSPOSE','BACKWARD','COLUMNWISE',m - k + i + ib - &
                    1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1),ldwork)

                 end if
                 ! apply h to rows 1:m-k+i+ib-1 of current block
                 call la_wung2l(m - k + i + ib - 1,ib,ib,a(1,n - k + i),lda,tau(i),work,iinfo &
                           )
                 ! set rows m-k+i+ib:m of current block to czero
                 do j = n - k + i,n - k + i + ib - 1
                    do l = m - k + i + ib,m
                       a(l,j) = czero
                    end do
                 end do
              end do
           end if
           work(1) = iws
           return
     end subroutine la_wungql
#endif

     pure subroutine la_cunm22(side,trans,m,n,n1,n2,q,ldq,c,ldc,work,lwork,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: m,n,n1,n2,ldq,ldc,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(in) :: q(ldq,*)
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,ldwork,len,lwkopt,nb,nq,nw
           ! Intrinsic Functions
           intrinsic :: cmplx,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q;
           ! nw is the minimum dimension of work.
           if (left) then
              nq = m
           else
              nq = n
           end if
           nw = nq
           if (n1 == 0 .or. n2 == 0) nw = 1
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (n1 < 0 .or. n1 + n2 /= nq) then
              info = -5
           else if (n2 < 0) then
              info = -6
           else if (ldq < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              lwkopt = m*n
              work(1) = cmplx(lwkopt,KIND=sp)
           end if
           if (info /= 0) then
              call la_xerbla('CUNM22',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           ! degenerate cases (n1 = 0 or n2 = 0) are handled using la_ctrmm.
           if (n1 == 0) then
              call la_ctrmm(side,'UPPER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           else if (n2 == 0) then
              call la_ctrmm(side,'LOWER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           end if
           ! compute the largest chunk size available from the workspace.
           nb = max(1,min(lwork,lwkopt)/nq)
           if (left) then
              if (notran) then
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q12.
                    call la_clacpy('ALL',n1,len,c(n2 + 1,i),ldc,work,ldwork)
                    call la_ctrmm('LEFT','LOWER','NO TRANSPOSE','NON-UNIT',n1,len,cone, &
                              q(1,n2 + 1),ldq,work,ldwork)
                    ! multiply top part of c by q11.
                    call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',n1,len,n2,cone,q,ldq, &
                              c(1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q21.
                    call la_clacpy('ALL',n2,len,c(1,i),ldc,work(n1 + 1),ldwork)

                    call la_ctrmm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',n2,len,cone, &
                              q(n1 + 1,1),ldq,work(n1 + 1),ldwork)
                    ! multiply bottom part of c by q22.
                    call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',n2,len,n1,cone,q(n1 + 1, &
                              n2 + 1),ldq,c(n2 + 1,i),ldc,cone,work(n1 + 1),ldwork)
                    ! copy everything back.
                    call la_clacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              else
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q21**h.
                    call la_clacpy('ALL',n2,len,c(n1 + 1,i),ldc,work,ldwork)
                    call la_ctrmm('LEFT','UPPER','CONJUGATE','NON-UNIT',n2,len,cone,q( &
                              n1 + 1,1),ldq,work,ldwork)
                    ! multiply top part of c by q11**h.
                    call la_cgemm('CONJUGATE','NO TRANSPOSE',n2,len,n1,cone,q,ldq,c( &
                              1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q12**h.
                    call la_clacpy('ALL',n1,len,c(1,i),ldc,work(n2 + 1),ldwork)

                    call la_ctrmm('LEFT','LOWER','CONJUGATE','NON-UNIT',n1,len,cone,q( &
                              1,n2 + 1),ldq,work(n2 + 1),ldwork)
                    ! multiply bottom part of c by q22**h.
                    call la_cgemm('CONJUGATE','NO TRANSPOSE',n1,len,n2,cone,q(n1 + 1,n2 + &
                              1),ldq,c(n1 + 1,i),ldc,cone,work(n2 + 1),ldwork)
                    ! copy everything back.
                    call la_clacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              end if
           else
              if (notran) then
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q21.
                    call la_clacpy('ALL',len,n2,c(i,n1 + 1),ldc,work,ldwork)
                    call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',len,n2,cone, &
                               q(n1 + 1,1),ldq,work,ldwork)
                    ! multiply left part of c by q11.
                    call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',len,n2,n1,cone,c(i,1) &
                              ,ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q12.
                    call la_clacpy('ALL',len,n1,c(i,1),ldc,work(1 + n2*ldwork), &
                              ldwork)
                    call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','NON-UNIT',len,n1,cone, &
                               q(1,n2 + 1),ldq,work(1 + n2*ldwork),ldwork)
                    ! multiply right part of c by q22.
                    call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',len,n1,n2,cone,c(i,n1 + &
                              1),ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n2*ldwork),ldwork)
                    ! copy everything back.
                    call la_clacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              else
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q12**h.
                    call la_clacpy('ALL',len,n1,c(i,n2 + 1),ldc,work,ldwork)
                    call la_ctrmm('RIGHT','LOWER','CONJUGATE','NON-UNIT',len,n1,cone,q( &
                               1,n2 + 1),ldq,work,ldwork)
                    ! multiply left part of c by q11**h.
                    call la_cgemm('NO TRANSPOSE','CONJUGATE',len,n1,n2,cone,c(i,1), &
                              ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q21**h.
                    call la_clacpy('ALL',len,n2,c(i,1),ldc,work(1 + n1*ldwork), &
                              ldwork)
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE','NON-UNIT',len,n2,cone,q( &
                               n1 + 1,1),ldq,work(1 + n1*ldwork),ldwork)
                    ! multiply right part of c by q22**h.
                    call la_cgemm('NO TRANSPOSE','CONJUGATE',len,n2,n1,cone,c(i,n2 + 1) &
                              ,ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n1*ldwork),ldwork)
                    ! copy everything back.
                    call la_clacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              end if
           end if
           work(1) = cmplx(lwkopt,KIND=sp)
           return
     end subroutine la_cunm22
     pure subroutine la_zunm22(side,trans,m,n,n1,n2,q,ldq,c,ldc,work,lwork,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: m,n,n1,n2,ldq,ldc,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(in) :: q(ldq,*)
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,ldwork,len,lwkopt,nb,nq,nw
           ! Intrinsic Functions
           intrinsic :: cmplx,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q;
           ! nw is the minimum dimension of work.
           if (left) then
              nq = m
           else
              nq = n
           end if
           nw = nq
           if (n1 == 0 .or. n2 == 0) nw = 1
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (n1 < 0 .or. n1 + n2 /= nq) then
              info = -5
           else if (n2 < 0) then
              info = -6
           else if (ldq < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              lwkopt = m*n
              work(1) = cmplx(lwkopt,KIND=dp)
           end if
           if (info /= 0) then
              call la_xerbla('ZUNM22',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           ! degenerate cases (n1 = 0 or n2 = 0) are handled using la_ztrmm.
           if (n1 == 0) then
              call la_ztrmm(side,'UPPER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           else if (n2 == 0) then
              call la_ztrmm(side,'LOWER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           end if
           ! compute the largest chunk size available from the workspace.
           nb = max(1,min(lwork,lwkopt)/nq)
           if (left) then
              if (notran) then
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q12.
                    call la_zlacpy('ALL',n1,len,c(n2 + 1,i),ldc,work,ldwork)
                    call la_ztrmm('LEFT','LOWER','NO TRANSPOSE','NON-UNIT',n1,len,cone, &
                              q(1,n2 + 1),ldq,work,ldwork)
                    ! multiply top part of c by q11.
                    call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',n1,len,n2,cone,q,ldq, &
                              c(1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q21.
                    call la_zlacpy('ALL',n2,len,c(1,i),ldc,work(n1 + 1),ldwork)

                    call la_ztrmm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',n2,len,cone, &
                              q(n1 + 1,1),ldq,work(n1 + 1),ldwork)
                    ! multiply bottom part of c by q22.
                    call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',n2,len,n1,cone,q(n1 + 1, &
                              n2 + 1),ldq,c(n2 + 1,i),ldc,cone,work(n1 + 1),ldwork)
                    ! copy everything back.
                    call la_zlacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              else
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q21**h.
                    call la_zlacpy('ALL',n2,len,c(n1 + 1,i),ldc,work,ldwork)
                    call la_ztrmm('LEFT','UPPER','CONJUGATE','NON-UNIT',n2,len,cone,q( &
                              n1 + 1,1),ldq,work,ldwork)
                    ! multiply top part of c by q11**h.
                    call la_zgemm('CONJUGATE','NO TRANSPOSE',n2,len,n1,cone,q,ldq,c( &
                              1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q12**h.
                    call la_zlacpy('ALL',n1,len,c(1,i),ldc,work(n2 + 1),ldwork)

                    call la_ztrmm('LEFT','LOWER','CONJUGATE','NON-UNIT',n1,len,cone,q( &
                              1,n2 + 1),ldq,work(n2 + 1),ldwork)
                    ! multiply bottom part of c by q22**h.
                    call la_zgemm('CONJUGATE','NO TRANSPOSE',n1,len,n2,cone,q(n1 + 1,n2 + &
                              1),ldq,c(n1 + 1,i),ldc,cone,work(n2 + 1),ldwork)
                    ! copy everything back.
                    call la_zlacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              end if
           else
              if (notran) then
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q21.
                    call la_zlacpy('ALL',len,n2,c(i,n1 + 1),ldc,work,ldwork)
                    call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',len,n2,cone, &
                               q(n1 + 1,1),ldq,work,ldwork)
                    ! multiply left part of c by q11.
                    call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',len,n2,n1,cone,c(i,1) &
                              ,ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q12.
                    call la_zlacpy('ALL',len,n1,c(i,1),ldc,work(1 + n2*ldwork), &
                              ldwork)
                    call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','NON-UNIT',len,n1,cone, &
                               q(1,n2 + 1),ldq,work(1 + n2*ldwork),ldwork)
                    ! multiply right part of c by q22.
                    call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',len,n1,n2,cone,c(i,n1 + &
                              1),ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n2*ldwork),ldwork)
                    ! copy everything back.
                    call la_zlacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              else
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q12**h.
                    call la_zlacpy('ALL',len,n1,c(i,n2 + 1),ldc,work,ldwork)
                    call la_ztrmm('RIGHT','LOWER','CONJUGATE','NON-UNIT',len,n1,cone,q( &
                               1,n2 + 1),ldq,work,ldwork)
                    ! multiply left part of c by q11**h.
                    call la_zgemm('NO TRANSPOSE','CONJUGATE',len,n1,n2,cone,c(i,1), &
                              ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q21**h.
                    call la_zlacpy('ALL',len,n2,c(i,1),ldc,work(1 + n1*ldwork), &
                              ldwork)
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE','NON-UNIT',len,n2,cone,q( &
                               n1 + 1,1),ldq,work(1 + n1*ldwork),ldwork)
                    ! multiply right part of c by q22**h.
                    call la_zgemm('NO TRANSPOSE','CONJUGATE',len,n2,n1,cone,c(i,n2 + 1) &
                              ,ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n1*ldwork),ldwork)
                    ! copy everything back.
                    call la_zlacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              end if
           end if
           work(1) = cmplx(lwkopt,KIND=dp)
           return
     end subroutine la_zunm22
#ifdef LA_WITH_XDP
     pure subroutine la_yunm22(side,trans,m,n,n1,n2,q,ldq,c,ldc,work,lwork,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: m,n,n1,n2,ldq,ldc,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(in) :: q(ldq,*)
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,ldwork,len,lwkopt,nb,nq,nw
           ! Intrinsic Functions
           intrinsic :: cmplx,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q;
           ! nw is the minimum dimension of work.
           if (left) then
              nq = m
           else
              nq = n
           end if
           nw = nq
           if (n1 == 0 .or. n2 == 0) nw = 1
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (n1 < 0 .or. n1 + n2 /= nq) then
              info = -5
           else if (n2 < 0) then
              info = -6
           else if (ldq < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              lwkopt = m*n
              work(1) = cmplx(lwkopt,KIND=xdp)
           end if
           if (info /= 0) then
              call la_xerbla('YUNM22',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           ! degenerate cases (n1 = 0 or n2 = 0) are handled using la_ytrmm.
           if (n1 == 0) then
              call la_ytrmm(side,'UPPER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           else if (n2 == 0) then
              call la_ytrmm(side,'LOWER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           end if
           ! compute the largest chunk size available from the workspace.
           nb = max(1,min(lwork,lwkopt)/nq)
           if (left) then
              if (notran) then
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q12.
                    call la_ylacpy('ALL',n1,len,c(n2 + 1,i),ldc,work,ldwork)
                    call la_ytrmm('LEFT','LOWER','NO TRANSPOSE','NON-UNIT',n1,len,cone, &
                              q(1,n2 + 1),ldq,work,ldwork)
                    ! multiply top part of c by q11.
                    call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',n1,len,n2,cone,q,ldq, &
                              c(1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q21.
                    call la_ylacpy('ALL',n2,len,c(1,i),ldc,work(n1 + 1),ldwork)

                    call la_ytrmm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',n2,len,cone, &
                              q(n1 + 1,1),ldq,work(n1 + 1),ldwork)
                    ! multiply bottom part of c by q22.
                    call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',n2,len,n1,cone,q(n1 + 1, &
                              n2 + 1),ldq,c(n2 + 1,i),ldc,cone,work(n1 + 1),ldwork)
                    ! copy everything back.
                    call la_ylacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              else
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q21**h.
                    call la_ylacpy('ALL',n2,len,c(n1 + 1,i),ldc,work,ldwork)
                    call la_ytrmm('LEFT','UPPER','CONJUGATE','NON-UNIT',n2,len,cone,q( &
                              n1 + 1,1),ldq,work,ldwork)
                    ! multiply top part of c by q11**h.
                    call la_ygemm('CONJUGATE','NO TRANSPOSE',n2,len,n1,cone,q,ldq,c( &
                              1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q12**h.
                    call la_ylacpy('ALL',n1,len,c(1,i),ldc,work(n2 + 1),ldwork)

                    call la_ytrmm('LEFT','LOWER','CONJUGATE','NON-UNIT',n1,len,cone,q( &
                              1,n2 + 1),ldq,work(n2 + 1),ldwork)
                    ! multiply bottom part of c by q22**h.
                    call la_ygemm('CONJUGATE','NO TRANSPOSE',n1,len,n2,cone,q(n1 + 1,n2 + &
                              1),ldq,c(n1 + 1,i),ldc,cone,work(n2 + 1),ldwork)
                    ! copy everything back.
                    call la_ylacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              end if
           else
              if (notran) then
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q21.
                    call la_ylacpy('ALL',len,n2,c(i,n1 + 1),ldc,work,ldwork)
                    call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',len,n2,cone, &
                               q(n1 + 1,1),ldq,work,ldwork)
                    ! multiply left part of c by q11.
                    call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',len,n2,n1,cone,c(i,1) &
                              ,ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q12.
                    call la_ylacpy('ALL',len,n1,c(i,1),ldc,work(1 + n2*ldwork), &
                              ldwork)
                    call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','NON-UNIT',len,n1,cone, &
                               q(1,n2 + 1),ldq,work(1 + n2*ldwork),ldwork)
                    ! multiply right part of c by q22.
                    call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',len,n1,n2,cone,c(i,n1 + &
                              1),ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n2*ldwork),ldwork)
                    ! copy everything back.
                    call la_ylacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              else
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q12**h.
                    call la_ylacpy('ALL',len,n1,c(i,n2 + 1),ldc,work,ldwork)
                    call la_ytrmm('RIGHT','LOWER','CONJUGATE','NON-UNIT',len,n1,cone,q( &
                               1,n2 + 1),ldq,work,ldwork)
                    ! multiply left part of c by q11**h.
                    call la_ygemm('NO TRANSPOSE','CONJUGATE',len,n1,n2,cone,c(i,1), &
                              ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q21**h.
                    call la_ylacpy('ALL',len,n2,c(i,1),ldc,work(1 + n1*ldwork), &
                              ldwork)
                    call la_ytrmm('RIGHT','UPPER','CONJUGATE','NON-UNIT',len,n2,cone,q( &
                               n1 + 1,1),ldq,work(1 + n1*ldwork),ldwork)
                    ! multiply right part of c by q22**h.
                    call la_ygemm('NO TRANSPOSE','CONJUGATE',len,n2,n1,cone,c(i,n2 + 1) &
                              ,ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n1*ldwork),ldwork)
                    ! copy everything back.
                    call la_ylacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              end if
           end if
           work(1) = cmplx(lwkopt,KIND=xdp)
           return
     end subroutine la_yunm22
#endif
#ifdef LA_WITH_QP
     pure subroutine la_wunm22(side,trans,m,n,n1,n2,q,ldq,c,ldc,work,lwork,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(in) :: m,n,n1,n2,ldq,ldc,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(in) :: q(ldq,*)
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,ldwork,len,lwkopt,nb,nq,nw
           ! Intrinsic Functions
           intrinsic :: cmplx,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q;
           ! nw is the minimum dimension of work.
           if (left) then
              nq = m
           else
              nq = n
           end if
           nw = nq
           if (n1 == 0 .or. n2 == 0) nw = 1
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'C')) &
                     then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (n1 < 0 .or. n1 + n2 /= nq) then
              info = -5
           else if (n2 < 0) then
              info = -6
           else if (ldq < max(1,nq)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              lwkopt = m*n
              work(1) = cmplx(lwkopt,KIND=qp)
           end if
           if (info /= 0) then
              call la_xerbla('WUNM22',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           ! degenerate cases (n1 = 0 or n2 = 0) are handled using la_wtrmm.
           if (n1 == 0) then
              call la_wtrmm(side,'UPPER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           else if (n2 == 0) then
              call la_wtrmm(side,'LOWER',trans,'NON-UNIT',m,n,cone,q,ldq,c,ldc)

              work(1) = cone
              return
           end if
           ! compute the largest chunk size available from the workspace.
           nb = max(1,min(lwork,lwkopt)/nq)
           if (left) then
              if (notran) then
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q12.
                    call la_wlacpy('ALL',n1,len,c(n2 + 1,i),ldc,work,ldwork)
                    call la_wtrmm('LEFT','LOWER','NO TRANSPOSE','NON-UNIT',n1,len,cone, &
                              q(1,n2 + 1),ldq,work,ldwork)
                    ! multiply top part of c by q11.
                    call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',n1,len,n2,cone,q,ldq, &
                              c(1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q21.
                    call la_wlacpy('ALL',n2,len,c(1,i),ldc,work(n1 + 1),ldwork)

                    call la_wtrmm('LEFT','UPPER','NO TRANSPOSE','NON-UNIT',n2,len,cone, &
                              q(n1 + 1,1),ldq,work(n1 + 1),ldwork)
                    ! multiply bottom part of c by q22.
                    call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',n2,len,n1,cone,q(n1 + 1, &
                              n2 + 1),ldq,c(n2 + 1,i),ldc,cone,work(n1 + 1),ldwork)
                    ! copy everything back.
                    call la_wlacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              else
                 do i = 1,n,nb
                    len = min(nb,n - i + 1)
                    ldwork = m
                    ! multiply bottom part of c by q21**h.
                    call la_wlacpy('ALL',n2,len,c(n1 + 1,i),ldc,work,ldwork)
                    call la_wtrmm('LEFT','UPPER','CONJUGATE','NON-UNIT',n2,len,cone,q( &
                              n1 + 1,1),ldq,work,ldwork)
                    ! multiply top part of c by q11**h.
                    call la_wgemm('CONJUGATE','NO TRANSPOSE',n2,len,n1,cone,q,ldq,c( &
                              1,i),ldc,cone,work,ldwork)
                    ! multiply top part of c by q12**h.
                    call la_wlacpy('ALL',n1,len,c(1,i),ldc,work(n2 + 1),ldwork)

                    call la_wtrmm('LEFT','LOWER','CONJUGATE','NON-UNIT',n1,len,cone,q( &
                              1,n2 + 1),ldq,work(n2 + 1),ldwork)
                    ! multiply bottom part of c by q22**h.
                    call la_wgemm('CONJUGATE','NO TRANSPOSE',n1,len,n2,cone,q(n1 + 1,n2 + &
                              1),ldq,c(n1 + 1,i),ldc,cone,work(n2 + 1),ldwork)
                    ! copy everything back.
                    call la_wlacpy('ALL',m,len,work,ldwork,c(1,i),ldc)
                 end do
              end if
           else
              if (notran) then
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q21.
                    call la_wlacpy('ALL',len,n2,c(i,n1 + 1),ldc,work,ldwork)
                    call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','NON-UNIT',len,n2,cone, &
                               q(n1 + 1,1),ldq,work,ldwork)
                    ! multiply left part of c by q11.
                    call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',len,n2,n1,cone,c(i,1) &
                              ,ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q12.
                    call la_wlacpy('ALL',len,n1,c(i,1),ldc,work(1 + n2*ldwork), &
                              ldwork)
                    call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','NON-UNIT',len,n1,cone, &
                               q(1,n2 + 1),ldq,work(1 + n2*ldwork),ldwork)
                    ! multiply right part of c by q22.
                    call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',len,n1,n2,cone,c(i,n1 + &
                              1),ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n2*ldwork),ldwork)
                    ! copy everything back.
                    call la_wlacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              else
                 do i = 1,m,nb
                    len = min(nb,m - i + 1)
                    ldwork = len
                    ! multiply right part of c by q12**h.
                    call la_wlacpy('ALL',len,n1,c(i,n2 + 1),ldc,work,ldwork)
                    call la_wtrmm('RIGHT','LOWER','CONJUGATE','NON-UNIT',len,n1,cone,q( &
                               1,n2 + 1),ldq,work,ldwork)
                    ! multiply left part of c by q11**h.
                    call la_wgemm('NO TRANSPOSE','CONJUGATE',len,n1,n2,cone,c(i,1), &
                              ldc,q,ldq,cone,work,ldwork)
                    ! multiply left part of c by q21**h.
                    call la_wlacpy('ALL',len,n2,c(i,1),ldc,work(1 + n1*ldwork), &
                              ldwork)
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE','NON-UNIT',len,n2,cone,q( &
                               n1 + 1,1),ldq,work(1 + n1*ldwork),ldwork)
                    ! multiply right part of c by q22**h.
                    call la_wgemm('NO TRANSPOSE','CONJUGATE',len,n2,n1,cone,c(i,n2 + 1) &
                              ,ldc,q(n1 + 1,n2 + 1),ldq,cone,work(1 + n1*ldwork),ldwork)
                    ! copy everything back.
                    call la_wlacpy('ALL',len,n,work,ldwork,c(i,1),ldc)
                 end do
              end if
           end if
           work(1) = cmplx(lwkopt,KIND=qp)
           return
     end subroutine la_wunm22
#endif

     !> CUNM2L: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by CGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_cunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           complex(sp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CUNM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = cone
              call la_clarf(side,mi,ni,a(1,i),1,taui,c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_cunm2l
     !> ZUNM2L: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by ZGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_zunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           complex(dp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZUNM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = cone
              call la_zlarf(side,mi,ni,a(1,i),1,taui,c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_zunm2l
#ifdef LA_WITH_XDP
     !> YUNM2L: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by YGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_yunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           complex(xdp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('YUNM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = cone
              call la_ylarf(side,mi,ni,a(1,i),1,taui,c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_yunm2l
#endif
#ifdef LA_WITH_QP
     !> WUNM2L: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by WGEQLF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_wunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,mi,ni,nq
           complex(qp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WUNM2L',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
           else
              mi = m
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(1:m-k+i,1:n)
                 mi = m - k + i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,1:n-k+i)
                 ni = n - k + i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              aii = a(nq - k + i,i)
              a(nq - k + i,i) = cone
              call la_wlarf(side,mi,ni,a(1,i),1,taui,c,ldc,work)
              a(nq - k + i,i) = aii
           end do
           return
     end subroutine la_wunm2l
#endif

     !> CUNML2: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by CGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_cunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           complex(sp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CUNML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = conjg(tau(i))
              else
                 taui = tau(i)
              end if
              if (i < nq) call la_clacgv(nq - i,a(i,i + 1),lda)
              aii = a(i,i)
              a(i,i) = cone
              call la_clarf(side,mi,ni,a(i,i),lda,taui,c(ic,jc),ldc,work)

              a(i,i) = aii
              if (i < nq) call la_clacgv(nq - i,a(i,i + 1),lda)
           end do
           return
     end subroutine la_cunml2
     !> ZUNML2: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by ZGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_zunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           complex(dp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZUNML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = conjg(tau(i))
              else
                 taui = tau(i)
              end if
              if (i < nq) call la_zlacgv(nq - i,a(i,i + 1),lda)
              aii = a(i,i)
              a(i,i) = cone
              call la_zlarf(side,mi,ni,a(i,i),lda,taui,c(ic,jc),ldc,work)

              a(i,i) = aii
              if (i < nq) call la_zlacgv(nq - i,a(i,i + 1),lda)
           end do
           return
     end subroutine la_zunml2
#ifdef LA_WITH_XDP
     !> YUNML2: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by YGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_yunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           complex(xdp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('YUNML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = conjg(tau(i))
              else
                 taui = tau(i)
              end if
              if (i < nq) call la_ylacgv(nq - i,a(i,i + 1),lda)
              aii = a(i,i)
              a(i,i) = cone
              call la_ylarf(side,mi,ni,a(i,i),lda,taui,c(ic,jc),ldc,work)

              a(i,i) = aii
              if (i < nq) call la_ylacgv(nq - i,a(i,i + 1),lda)
           end do
           return
     end subroutine la_yunml2
#endif
#ifdef LA_WITH_QP
     !> WUNML2: overwrites the general complex m-by-n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by WGELQF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_wunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,jc,mi,ni,nq
           complex(qp) :: aii,taui
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           ! nq is the order of q
           if (left) then
              nq = m
           else
              nq = n
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WUNML2',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. notran .or. .not. left .and. .not. notran)) then
              i1 = 1
              i2 = k
              i3 = 1
           else
              i1 = k
              i2 = 1
              i3 = -1
           end if
           if (left) then
              ni = n
              jc = 1
           else
              mi = m
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**h is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**h is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**h
              if (notran) then
                 taui = conjg(tau(i))
              else
                 taui = tau(i)
              end if
              if (i < nq) call la_wlacgv(nq - i,a(i,i + 1),lda)
              aii = a(i,i)
              a(i,i) = cone
              call la_wlarf(side,mi,ni,a(i,i),lda,taui,c(ic,jc),ldc,work)

              a(i,i) = aii
              if (i < nq) call la_wlacgv(nq - i,a(i,i + 1),lda)
           end do
           return
     end subroutine la_wunml2
#endif

     !> CUNMLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by CGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_cunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0 .or. k == 0) then
                 lwkopt = 1
              else
              nb = min(nbmax,la_ilaenv(1,'CUNMLQ',side//trans,m,n,k,-1))

              lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              return
           end if
           ! determine the block size
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'CUNMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_cunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_clarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**h is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**h
                 call la_clarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_cunmlq
     !> ZUNMLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by ZGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_zunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'ZUNMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'ZUNMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_zunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_zlarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**h is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**h
                 call la_zlarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_zunmlq
#ifdef LA_WITH_XDP
     !> YUNMLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by YGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_yunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'YUNMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YUNMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'YUNMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_yunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_ylarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**h is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**h
                 call la_ylarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_yunmlq
#endif
#ifdef LA_WITH_QP
     !> WUNMLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k)**H . . . H(2)**H H(1)**H
     !> as returned by WGELQF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_wunmlq(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           character :: transt
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,jc,ldwork,lwkopt,mi,nb,nbmin, &
                     ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,k)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              nb = min(nbmax,la_ilaenv(1,'WUNMLQ',side//trans,m,n,k,-1))
              lwkopt = nw*nb + tsize
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNMLQ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'WUNMLQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_wunml2(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
                 jc = 1
              else
                 mi = m
                 ic = 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i) h(i+1) . . . h(i+ib-1)
                 call la_wlarft('FORWARD','ROWWISE',nq - i + 1,ib,a(i,i),lda,tau(i), &
                           work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(i:m,1:n)
                    mi = m - i + 1
                    ic = i
                 else
                    ! h or h**h is applied to c(1:m,i:n)
                    ni = n - i + 1
                    jc = i
                 end if
                 ! apply h or h**h
                 call la_wlarfb(side,transt,'FORWARD','ROWWISE',mi,ni,ib,a(i,i), &
                           lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_wunmlq
#endif

     !> CUNMQL: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by CGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_cunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'CUNMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           ! determine the block size
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'CUNMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_cunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_clarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**h is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**h
                 call la_clarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_cunmql
     !> ZUNMQL: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by ZGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_zunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'ZUNMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'ZUNMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_zunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_zlarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**h is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**h
                 call la_zlarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_zunmql
#ifdef LA_WITH_XDP
     !> YUNMQL: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by YGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_yunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'YUNMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YUNMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'YUNMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_yunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_ylarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**h is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**h
                 call la_ylarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_yunmql
#endif
#ifdef LA_WITH_QP
     !> WUNMQL: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(k) . . . H(2) H(1)
     !> as returned by WGEQLF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_wunmql(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: nbmax = 64
           integer(ilp),parameter :: ldt = nbmax + 1
           integer(ilp),parameter :: tsize = ldt*nbmax

           ! Local Scalars
           logical(lk) :: left,lquery,notran
           integer(ilp) :: i,i1,i2,i3,ib,iinfo,iwt,ldwork,lwkopt,mi,nb,nbmin,ni,nq, &
                     nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > nq) then
              info = -5
           else if (lda < max(1,nq)) then
              info = -7
           else if (ldc < max(1,m)) then
              info = -10
           else if (lwork < nw .and. .not. lquery) then
              info = -12
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'WUNMQL',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNMQL',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'WUNMQL',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_wunm2l(side,trans,m,n,k,a,lda,tau,c,ldc,work,iinfo)
           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. notran) .or. (.not. left .and. .not. notran)) then
                 i1 = 1
                 i2 = k
                 i3 = nb
              else
                 i1 = ((k - 1)/nb)*nb + 1
                 i2 = 1
                 i3 = -nb
              end if
              if (left) then
                 ni = n
              else
                 mi = m
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_wlarft('BACKWARD','COLUMNWISE',nq - k + i + ib - 1,ib,a(1,i),lda, &
                           tau(i),work(iwt),ldt)
                 if (left) then
                    ! h or h**h is applied to c(1:m-k+i+ib-1,1:n)
                    mi = m - k + i + ib - 1
                 else
                    ! h or h**h is applied to c(1:m,1:n-k+i+ib-1)
                    ni = n - k + i + ib - 1
                 end if
                 ! apply h or h**h
                 call la_wlarfb(side,trans,'BACKWARD','COLUMNWISE',mi,ni,ib,a(1,i), &
                           lda,work(iwt),ldt,c,ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_wunmql
#endif

     !> CGELQ2: computes an LQ factorization of a complex m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_cgelq2(m,n,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_clacgv(n - i + 1,a(i,i),lda)
              alpha = a(i,i)
              call la_clarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 a(i,i) = cone
                 call la_clarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
              end if
              a(i,i) = alpha
              call la_clacgv(n - i + 1,a(i,i),lda)
           end do
           return
     end subroutine la_cgelq2
     !> ZGELQ2: computes an LQ factorization of a complex m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_zgelq2(m,n,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_zlacgv(n - i + 1,a(i,i),lda)
              alpha = a(i,i)
              call la_zlarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 a(i,i) = cone
                 call la_zlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
              end if
              a(i,i) = alpha
              call la_zlacgv(n - i + 1,a(i,i),lda)
           end do
           return
     end subroutine la_zgelq2
#ifdef LA_WITH_XDP
     !> YGELQ2: computes an LQ factorization of a complex m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_ygelq2(m,n,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('YGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_ylacgv(n - i + 1,a(i,i),lda)
              alpha = a(i,i)
              call la_ylarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 a(i,i) = cone
                 call la_ylarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
              end if
              a(i,i) = alpha
              call la_ylacgv(n - i + 1,a(i,i),lda)
           end do
           return
     end subroutine la_ygelq2
#endif
#ifdef LA_WITH_QP
     !> WGELQ2: computes an LQ factorization of a complex m-by-n matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a n-by-n orthogonal matrix;
     !> L is a lower-triangular m-by-m matrix;
     !> 0 is a m-by-(n-m) zero matrix, if m < n.

     pure subroutine la_wgelq2(m,n,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WGELQ2',-info)
              return
           end if
           k = min(m,n)
           do i = 1,k
              ! generate elementary reflector h(i) to annihilate a(i,i+1:n)
              call la_wlacgv(n - i + 1,a(i,i),lda)
              alpha = a(i,i)
              call la_wlarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,tau(i))
              if (i < m) then
                 ! apply h(i) to a(i+1:m,i:n) from the right
                 a(i,i) = cone
                 call la_wlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,tau(i),a(i + 1,i), &
                           lda,work)
              end if
              a(i,i) = alpha
              call la_wlacgv(n - i + 1,a(i,i),lda)
           end do
           return
     end subroutine la_wgelq2
#endif

     !> CGELQF: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_cgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'CGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'CGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'CGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_cgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_clarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_clarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_cgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_cgelqf
     !> ZGELQF: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_zgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'ZGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'ZGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'ZGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_zgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_zlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_zlarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_zgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_zgelqf
#ifdef LA_WITH_XDP
     !> YGELQF: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_ygelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'YGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'YGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'YGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_ygelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_ylarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_ylarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_ygelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_ygelqf
#endif
#ifdef LA_WITH_QP
     !> WGELQF: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_wgelqf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ldwork,lwkopt,nb,nbmin,nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           nb = la_ilaenv(1,'WGELQF',' ',m,n,-1,-1)
           lwkopt = m*nb
           work(1) = lwkopt
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m) .and. .not. lquery) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WGELQF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           nx = 0
           iws = m
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'WGELQF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'WGELQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially
              do i = 1,k - nx,nb
                 ib = min(k - i + 1,nb)
                 ! compute the lq factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_wgelq2(ib,n - i + 1,a(i,i),lda,tau(i),work,iinfo)
                 if (i + ib <= m) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i) h(i+1) . . . h(i+ib-1)
                    call la_wlarft('FORWARD','ROWWISE',n - i + 1,ib,a(i,i),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(i+ib:m,i:n) from the right
                    call la_wlarfb('RIGHT','NO TRANSPOSE','FORWARD','ROWWISE',m - i - ib + 1,n - &
                    i + 1,ib,a(i,i),lda,work,ldwork,a(i + ib,i),lda,work(ib + 1),ldwork)

                 end if
              end do
           else
              i = 1
           end if
           ! use unblocked code to factor the last or only block.
           if (i <= k) call la_wgelq2(m - i + 1,n - i + 1,a(i,i),lda,tau(i),work,iinfo)

           work(1) = iws
           return
     end subroutine la_wgelqf
#endif

     !> CGELQT3: recursively computes a LQ factorization of a complex M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_cgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_clarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
              t(1,1) = conjg(t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_cgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) =  a(j1:m,1:n) q1^h [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_ctrmm('R','U','C','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              call la_cgemm('N','C',m2,m1,n - m1,cone,a(i1,i1),lda,a(1,i1),lda, &
                        cone,t(i1,1),ldt)
              call la_ctrmm('R','U','N','N',m2,m1,cone,t,ldt,t(i1,1),ldt)

              call la_cgemm('N','N',m2,n - m1,m1,-cone,t(i1,1),ldt,a(1,i1),lda, &
                        cone,a(i1,i1),lda)
              call la_ctrmm('R','U','N','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = czero
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_cgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_ctrmm('R','U','C','U',m1,m2,cone,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_cgemm('N','C',m1,m2,n - m,cone,a(1,j1),lda,a(i1,j1),lda, &
                        cone,t(1,i1),ldt)
              call la_ctrmm('L','U','N','N',m1,m2,-cone,t,ldt,t(1,i1),ldt)

              call la_ctrmm('R','U','N','N',m1,m2,cone,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_cgelqt3
     !> ZGELQT3: recursively computes a LQ factorization of a complex M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_zgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_zlarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
              t(1,1) = conjg(t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_zgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) =  a(j1:m,1:n) q1^h [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_ztrmm('R','U','C','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              call la_zgemm('N','C',m2,m1,n - m1,cone,a(i1,i1),lda,a(1,i1),lda, &
                        cone,t(i1,1),ldt)
              call la_ztrmm('R','U','N','N',m2,m1,cone,t,ldt,t(i1,1),ldt)

              call la_zgemm('N','N',m2,n - m1,m1,-cone,t(i1,1),ldt,a(1,i1),lda, &
                        cone,a(i1,i1),lda)
              call la_ztrmm('R','U','N','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = czero
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_zgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_ztrmm('R','U','C','U',m1,m2,cone,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_zgemm('N','C',m1,m2,n - m,cone,a(1,j1),lda,a(i1,j1),lda, &
                        cone,t(1,i1),ldt)
              call la_ztrmm('L','U','N','N',m1,m2,-cone,t,ldt,t(1,i1),ldt)

              call la_ztrmm('R','U','N','N',m1,m2,cone,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_zgelqt3
#ifdef LA_WITH_XDP
     !> YGELQT3: recursively computes a LQ factorization of a complex M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_ygelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('YGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_ylarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
              t(1,1) = conjg(t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_ygelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) =  a(j1:m,1:n) q1^h [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_ytrmm('R','U','C','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              call la_ygemm('N','C',m2,m1,n - m1,cone,a(i1,i1),lda,a(1,i1),lda, &
                        cone,t(i1,1),ldt)
              call la_ytrmm('R','U','N','N',m2,m1,cone,t,ldt,t(i1,1),ldt)

              call la_ygemm('N','N',m2,n - m1,m1,-cone,t(i1,1),ldt,a(1,i1),lda, &
                        cone,a(i1,i1),lda)
              call la_ytrmm('R','U','N','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = czero
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_ygelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_ytrmm('R','U','C','U',m1,m2,cone,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_ygemm('N','C',m1,m2,n - m,cone,a(1,j1),lda,a(i1,j1),lda, &
                        cone,t(1,i1),ldt)
              call la_ytrmm('L','U','N','N',m1,m2,-cone,t,ldt,t(1,i1),ldt)

              call la_ytrmm('R','U','N','N',m1,m2,cone,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_ygelqt3
#endif
#ifdef LA_WITH_QP
     !> WGELQT3: recursively computes a LQ factorization of a complex M-by-N
     !> matrix A, using the compact WY representation of Q.
     !> Based on the algorithm of Elmroth and Gustavson,
     !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.

     pure recursive subroutine la_wgelqt3(m,n,a,lda,t,ldt,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,ldt
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: t(ldt,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,j,j1,m1,m2,iinfo
           ! Executable Statements
           info = 0
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (ldt < max(1,m)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WGELQT3',-info)
              return
           end if
           if (m == 1) then
              ! compute householder transform when m=1
              call la_wlarfg(n,a(1,1),a(1,min(2,n)),lda,t(1,1))
              t(1,1) = conjg(t(1,1))
           else
              ! otherwise, split a into blocks...
              m1 = m/2
              m2 = m - m1
              i1 = min(m1 + 1,m)
              j1 = min(m + 1,n)
              ! compute a(1:m1,1:n) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
              call la_wgelqt3(m1,n,a,lda,t,ldt,iinfo)
              ! compute a(j1:m,1:n) =  a(j1:m,1:n) q1^h [workspace: t(1:n1,j1:n)]
              do i = 1,m2
                 do j = 1,m1
                    t(i + m1,j) = a(i + m1,j)
                 end do
              end do
              call la_wtrmm('R','U','C','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              call la_wgemm('N','C',m2,m1,n - m1,cone,a(i1,i1),lda,a(1,i1),lda, &
                        cone,t(i1,1),ldt)
              call la_wtrmm('R','U','N','N',m2,m1,cone,t,ldt,t(i1,1),ldt)

              call la_wgemm('N','N',m2,n - m1,m1,-cone,t(i1,1),ldt,a(1,i1),lda, &
                        cone,a(i1,i1),lda)
              call la_wtrmm('R','U','N','U',m2,m1,cone,a,lda,t(i1,1),ldt)

              do i = 1,m2
                 do j = 1,m1
                    a(i + m1,j) = a(i + m1,j) - t(i + m1,j)
                    t(i + m1,j) = czero
                 end do
              end do
              ! compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
              call la_wgelqt3(m2,n - m1,a(i1,i1),lda,t(i1,i1),ldt,iinfo)
              ! compute t3 = t(j1:n1,1:n) = -t1 y1^h y2 t2
              do i = 1,m2
                 do j = 1,m1
                    t(j,i + m1) = (a(j,i + m1))
                 end do
              end do
              call la_wtrmm('R','U','C','U',m1,m2,cone,a(i1,i1),lda,t(1,i1), &
                        ldt)
              call la_wgemm('N','C',m1,m2,n - m,cone,a(1,j1),lda,a(i1,j1),lda, &
                        cone,t(1,i1),ldt)
              call la_wtrmm('L','U','N','N',m1,m2,-cone,t,ldt,t(1,i1),ldt)

              call la_wtrmm('R','U','N','N',m1,m2,cone,t(i1,i1),ldt,t(1,i1), &
                        ldt)
              ! y = (y1,y2); l = [ l1            0  ];  t = [t1 t3]
                               ! [ a(1:n1,j1:n)  l2 ]       [ 0 t2]
           end if
           return
     end subroutine la_wgelqt3
#endif

     !> CGEMLQT: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'C':   Q**H C            C Q**H
     !> where Q is a complex unitary matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**H
     !> generated using the compact WY representation as returned by CGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_cgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           complex(sp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_clarfb('L','C','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_clarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_clarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_clarfb('R','C','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_cgemlqt
     !> ZGEMLQT: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'C':   Q**H C            C Q**H
     !> where Q is a complex unitary matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**H
     !> generated using the compact WY representation as returned by ZGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_zgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           complex(dp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_zlarfb('L','C','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_zlarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_zlarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_zlarfb('R','C','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_zgemlqt
#ifdef LA_WITH_XDP
     !> YGEMLQT: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'C':   Q**H C            C Q**H
     !> where Q is a complex unitary matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**H
     !> generated using the compact WY representation as returned by YGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_ygemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           complex(xdp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('YGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_ylarfb('L','C','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_ylarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_ylarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_ylarfb('R','C','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_ygemlqt
#endif
#ifdef LA_WITH_QP
     !> WGEMLQT: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q C            C Q
     !> TRANS = 'C':   Q**H C            C Q**H
     !> where Q is a complex unitary matrix defined as the product of K
     !> elementary reflectors:
     !> Q = H(1) H(2) . . . H(K) = I - V T V**H
     !> generated using the compact WY representation as returned by WGELQT.
     !> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.

     pure subroutine la_wgemlqt(side,trans,m,n,k,mb,v,ldv,t,ldt,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,ldc,m,n,mb,ldt
           ! Array Arguments
           complex(qp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,ldwork,kf,q
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldwork = max(1,n)
              q = m
           else if (right) then
              ldwork = max(1,m)
              q = n
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0 .or. k > q) then
              info = -5
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -6
           else if (ldv < max(1,k)) then
               info = -8
           else if (ldt < mb) then
              info = -10
           else if (ldc < max(1,m)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WGEMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_wlarfb('L','C','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 call la_wlarfb('R','N','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_wlarfb('L','N','F','R',m - i + 1,n,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(i,1),ldc,work,ldwork)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 call la_wlarfb('R','C','F','R',m,n - i + 1,ib,v(i,i),ldv,t(1,i), &
                           ldt,c(1,i),ldc,work,ldwork)
              end do
           end if
           return
     end subroutine la_wgemlqt
#endif

     !> CGEQL2: computes a QL factorization of a complex m by n matrix A:
     !> A = Q * L.

     pure subroutine la_cgeql2(m,n,a,lda,tau,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              alpha = a(m - k + i,n - k + i)
              call la_clarfg(m - k + i,alpha,a(1,n - k + i),1,tau(i))
              ! apply h(i)**h to a(1:m-k+i,1:n-k+i-1) from the left
              a(m - k + i,n - k + i) = cone
              call la_clarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,conjg(tau(i)),a, &
                        lda,work)
              a(m - k + i,n - k + i) = alpha
           end do
           return
     end subroutine la_cgeql2
     !> ZGEQL2: computes a QL factorization of a complex m by n matrix A:
     !> A = Q * L.

     pure subroutine la_zgeql2(m,n,a,lda,tau,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              alpha = a(m - k + i,n - k + i)
              call la_zlarfg(m - k + i,alpha,a(1,n - k + i),1,tau(i))
              ! apply h(i)**h to a(1:m-k+i,1:n-k+i-1) from the left
              a(m - k + i,n - k + i) = cone
              call la_zlarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,conjg(tau(i)),a, &
                        lda,work)
              a(m - k + i,n - k + i) = alpha
           end do
           return
     end subroutine la_zgeql2
#ifdef LA_WITH_XDP
     !> YGEQL2: computes a QL factorization of a complex m by n matrix A:
     !> A = Q * L.

     pure subroutine la_ygeql2(m,n,a,lda,tau,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('YGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              alpha = a(m - k + i,n - k + i)
              call la_ylarfg(m - k + i,alpha,a(1,n - k + i),1,tau(i))
              ! apply h(i)**h to a(1:m-k+i,1:n-k+i-1) from the left
              a(m - k + i,n - k + i) = cone
              call la_ylarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,conjg(tau(i)),a, &
                        lda,work)
              a(m - k + i,n - k + i) = alpha
           end do
           return
     end subroutine la_ygeql2
#endif
#ifdef LA_WITH_QP
     !> WGEQL2: computes a QL factorization of a complex m by n matrix A:
     !> A = Q * L.

     pure subroutine la_wgeql2(m,n,a,lda,tau,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,k
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WGEQL2',-info)
              return
           end if
           k = min(m,n)
           do i = k,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! a(1:m-k+i-1,n-k+i)
              alpha = a(m - k + i,n - k + i)
              call la_wlarfg(m - k + i,alpha,a(1,n - k + i),1,tau(i))
              ! apply h(i)**h to a(1:m-k+i,1:n-k+i-1) from the left
              a(m - k + i,n - k + i) = cone
              call la_wlarf('LEFT',m - k + i,n - k + i - 1,a(1,n - k + i),1,conjg(tau(i)),a, &
                        lda,work)
              a(m - k + i,n - k + i) = alpha
           end do
           return
     end subroutine la_wgeql2
#endif

     !> CGEQLF: computes a QL factorization of a complex M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_cgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'CGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'CGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'CGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_cgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_clarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_clarfb('LEFT','CONJUGATE TRANSPOSE','BACKWARD','COLUMNWISE',m - &
                    k + i + ib - 1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1), &
                              ldwork)
                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_cgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_cgeqlf
     !> ZGEQLF: computes a QL factorization of a complex M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_zgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'ZGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'ZGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'ZGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_zgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_zlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_zlarfb('LEFT','CONJUGATE TRANSPOSE','BACKWARD','COLUMNWISE',m - &
                    k + i + ib - 1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1), &
                              ldwork)
                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_zgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_zgeqlf
#ifdef LA_WITH_XDP
     !> YGEQLF: computes a QL factorization of a complex M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_ygeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'YGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'YGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'YGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_ygeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_ylarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_ylarfb('LEFT','CONJUGATE TRANSPOSE','BACKWARD','COLUMNWISE',m - &
                    k + i + ib - 1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1), &
                              ldwork)
                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_ygeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_ygeqlf
#endif
#ifdef LA_WITH_QP
     !> WGEQLF: computes a QL factorization of a complex M-by-N matrix A:
     !> A = Q * L.

     pure subroutine la_wgeqlf(m,n,a,lda,tau,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ib,iinfo,iws,k,ki,kk,ldwork,lwkopt,mu,nb,nbmin,nu, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              k = min(m,n)
              if (k == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'WGEQLF',' ',m,n,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
              if (lwork < max(1,n) .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGEQLF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (k == 0) then
              return
           end if
           nbmin = 2
           nx = 1
           iws = n
           if (nb > 1 .and. nb < k) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'WGEQLF',' ',m,n,-1,-1))
              if (nx < k) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = n
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'WGEQLF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < k .and. nx < k) then
              ! use blocked code initially.
              ! the last kk columns are handled by the block method.
              ki = ((k - nx - 1)/nb)*nb
              kk = min(k,ki + nb)
              do i = k - kk + ki + 1,k - kk + 1,-nb
                 ib = min(k - i + 1,nb)
                 ! compute the ql factorization of the current block
                 ! a(1:m-k+i+ib-1,n-k+i:n-k+i+ib-1)
                 call la_wgeql2(m - k + i + ib - 1,ib,a(1,n - k + i),lda,tau(i),work,iinfo)

                 if (n - k + i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_wlarft('BACKWARD','COLUMNWISE',m - k + i + ib - 1,ib,a(1,n - k + i), &
                              lda,tau(i),work,ldwork)
                    ! apply h**h to a(1:m-k+i+ib-1,1:n-k+i-1) from the left
                    call la_wlarfb('LEFT','CONJUGATE TRANSPOSE','BACKWARD','COLUMNWISE',m - &
                    k + i + ib - 1,n - k + i - 1,ib,a(1,n - k + i),lda,work,ldwork,a,lda,work(ib + 1), &
                              ldwork)
                 end if
              end do
              mu = m - k + i + nb - 1
              nu = n - k + i + nb - 1
           else
              mu = m
              nu = n
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0 .and. nu > 0) call la_wgeql2(mu,nu,a,lda,tau,work,iinfo)
           work(1) = iws
           return
     end subroutine la_wgeqlf
#endif

     !> CTPLQT: computes a blocked LQ factorization of a complex
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_ctplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_ctplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_ctprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_ctplqt
     !> ZTPLQT: computes a blocked LQ factorization of a complex
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_ztplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_ztplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_ztprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_ztplqt
#ifdef LA_WITH_XDP
     !> YTPLQT: computes a blocked LQ factorization of a complex
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_ytplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('YTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_ytplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_ytprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_ytplqt
#endif
#ifdef LA_WITH_QP
     !> WTPLQT: computes a blocked LQ factorization of a complex
     !> "triangular-pentagonal" matrix C, which is composed of a
     !> triangular block A and pentagonal block B, using the compact
     !> WY representation for Q.

     pure subroutine la_wtplqt(m,n,l,mb,a,lda,b,ldb,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldt,n,m,l,mb
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,lb,nb,iinfo
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (l < 0 .or. (l > min(m,n) .and. min(m,n) >= 0)) then
              info = -3
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (ldb < max(1,m)) then
              info = -8
           else if (ldt < mb) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WTPLQT',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           do i = 1,m,mb
           ! compute the qr factorization of the current block
              ib = min(m - i + 1,mb)
              nb = min(n - l + i + ib - 1,n)
              if (i >= l) then
                 lb = 0
              else
                 lb = nb - n + l - i + 1
              end if
              call la_wtplqt2(ib,nb,lb,a(i,i),lda,b(i,1),ldb,t(1,i),ldt,iinfo)

           ! update by applying h**t to b(i+ib:m,:) from the right
              if (i + ib <= m) then
                 call la_wtprfb('R','N','F','R',m - i - ib + 1,nb,ib,lb,b(i,1),ldb,t( &
                           1,i),ldt,a(i + ib,i),lda,b(i + ib,1),ldb,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_wtplqt
#endif

     !> CTPMLQT: applies a complex unitary matrix Q obtained from a
     !> "triangular-pentagonal" complex block reflector H to a general
     !> complex matrix C, which consists of two blocks A and B.

     pure subroutine la_ctpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           complex(sp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('CTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ctprfb('L','C','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ctprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ctprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ctprfb('R','C','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_ctpmlqt
     !> ZTPMLQT: applies a complex unitary matrix Q obtained from a
     !> "triangular-pentagonal" complex block reflector H to a general
     !> complex matrix C, which consists of two blocks A and B.

     pure subroutine la_ztpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           complex(dp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('ZTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ztprfb('L','C','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ztprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ztprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ztprfb('R','C','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_ztpmlqt
#ifdef LA_WITH_XDP
     !> YTPMLQT: applies a complex unitary matrix Q obtained from a
     !> "triangular-pentagonal" complex block reflector H to a general
     !> complex matrix C, which consists of two blocks A and B.

     pure subroutine la_ytpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           complex(xdp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('YTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ytprfb('L','C','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ytprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_ytprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_ytprfb('R','C','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_ytpmlqt
#endif
#ifdef LA_WITH_QP
     !> WTPMLQT: applies a complex unitary matrix Q obtained from a
     !> "triangular-pentagonal" complex block reflector H to a general
     !> complex matrix C, which consists of two blocks A and B.

     pure subroutine la_wtpmlqt(side,trans,m,n,k,l,mb,v,ldv,t,ldt,a,lda,b,ldb, &
               work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldv,lda,ldb,m,n,l,mb,ldt
           ! Array Arguments
           complex(qp),intent(in) :: v(ldv,*),t(ldt,*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran
           integer(ilp) :: i,ib,nb,lb,kf,ldaq
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! Test The Input Arguments
           info = 0
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           tran = la_lsame(trans,'C')
           notran = la_lsame(trans,'N')
           if (left) then
              ldaq = max(1,k)
           else if (right) then
              ldaq = max(1,m)
           end if
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (m < 0) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (l < 0 .or. l > k) then
              info = -6
           else if (mb < 1 .or. (mb > k .and. k > 0)) then
              info = -7
           else if (ldv < k) then
              info = -9
           else if (ldt < mb) then
              info = -11
           else if (lda < ldaq) then
              info = -13
           else if (ldb < max(1,m)) then
              info = -15
           end if
           if (info /= 0) then
              call la_xerbla('WTPMLQT',-info)
              return
           end if
           ! Quick Return If Possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if (left .and. notran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_wtprfb('L','C','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. tran) then
              do i = 1,k,mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_wtprfb('R','N','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           else if (left .and. tran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(m - l + i + ib - 1,m)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = 0
                 end if
                 call la_wtprfb('L','N','F','R',nb,n,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(i,1),lda,b,ldb,work,ib)
              end do
           else if (right .and. notran) then
              kf = ((k - 1)/mb)*mb + 1
              do i = kf,1,-mb
                 ib = min(mb,k - i + 1)
                 nb = min(n - l + i + ib - 1,n)
                 if (i >= l) then
                    lb = 0
                 else
                    lb = nb - n + l - i + 1
                 end if
                 call la_wtprfb('R','C','F','R',m,nb,ib,lb,v(i,1),ldv,t(1,i), &
                           ldt,a(1,i),lda,b,ldb,work,m)
              end do
           end if
           return
     end subroutine la_wtpmlqt
#endif

     !> CGELQT: computes a blocked LQ factorization of a complex M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_cgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_cgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_clarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_cgelqt
     !> ZGELQT: computes a blocked LQ factorization of a complex M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_zgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_zgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_zlarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_zgelqt
#ifdef LA_WITH_XDP
     !> YGELQT: computes a blocked LQ factorization of a complex M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_ygelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_ygelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_ylarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_ygelqt
#endif
#ifdef LA_WITH_QP
     !> WGELQT: computes a blocked LQ factorization of a complex M-by-N matrix A
     !> using the compact WY representation of Q.

     pure subroutine la_wgelqt(m,n,mb,a,lda,t,ldt,work,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldt,m,n,mb
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: t(ldt,*),work(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ib,iinfo,k
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (mb < 1 .or. (mb > min(m,n) .and. min(m,n) > 0)) then
              info = -3
           else if (lda < max(1,m)) then
              info = -5
           else if (ldt < mb) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WGELQT',-info)
              return
           end if
           ! quick return if possible
           k = min(m,n)
           if (k == 0) return
           ! blocked loop of length k
           do i = 1,k,mb
              ib = min(k - i + 1,mb)
           ! compute the lq factorization of the current block a(i:m,i:i+ib-1)
              call la_wgelqt3(ib,n - i + 1,a(i,i),lda,t(1,i),ldt,iinfo)
              if (i + ib <= m) then
           ! update by applying h**t to a(i:m,i+ib:n) from the right
              call la_wlarfb('R','N','F','R',m - i - ib + 1,n - i + 1,ib,a(i,i),lda,t(1,i &
                        ),ldt,a(i + ib,i),lda,work,m - i - ib + 1)
              end if
           end do
           return
     end subroutine la_wgelqt
#endif

     !> CLAMSWLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (CLASWLQ)

     pure subroutine la_clamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),t(ldt,*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,lw,ctr
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('CLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_cgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_ctpmlqt('L','C',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_ctpmlqt('L','C',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                           1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_cgemlqt('L','C',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_cgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_ctpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c( &
                         1,1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_ctpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_ctpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               ctr = ctr - 1
                   call la_ctpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                              c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_cgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ii = n - kk + 1
              ctr = 1
              call la_cgemlqt('R','C',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_ctpmlqt('R','C',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_ctpmlqt('R','C',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c(1,1), &
                          ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_clamswlq
     !> ZLAMSWLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (ZLASWLQ)

     pure subroutine la_zlamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),t(ldt,*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,lw,ctr
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('ZLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_zgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_ztpmlqt('L','C',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_ztpmlqt('L','C',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                           1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_zgemlqt('L','C',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_zgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_ztpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_ztpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_ztpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1), &
                           ldt,c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               ctr = ctr - 1
               call la_ztpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1), &
                         ldt,c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_zgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ii = n - kk + 1
              call la_zgemlqt('R','C',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              ctr = 1
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_ztpmlqt('R','C',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_ztpmlqt('R','C',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                         1,1),ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_zlamswlq
#ifdef LA_WITH_XDP
     !> YLAMSWLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (YLASWLQ)

     pure subroutine la_ylamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),t(ldt,*)
           complex(xdp),intent(out) :: work(*)
           complex(xdp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,lw,ctr
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('YLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_ygemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_ytpmlqt('L','C',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_ytpmlqt('L','C',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                           1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_ygemlqt('L','C',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_ygemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_ytpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_ytpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_ytpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1), &
                           ldt,c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               ctr = ctr - 1
               call la_ytpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1), &
                         ldt,c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_ygemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ii = n - kk + 1
              call la_ygemlqt('R','C',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              ctr = 1
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_ytpmlqt('R','C',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_ytpmlqt('R','C',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                         1,1),ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_ylamswlq
#endif
#ifdef LA_WITH_QP
     !> WLAMSWLQ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of blocked
     !> elementary reflectors computed by short wide LQ
     !> factorization (WLASWLQ)

     pure subroutine la_wlamswlq(side,trans,m,n,k,mb,nb,a,lda,t,ldt,c,ldc,work, &
               lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,mb,nb,ldt,lwork,ldc
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),t(ldt,*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: c(ldc,*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: i,ii,kk,lw,ctr
           ! External Subroutines
           ! Executable Statements
           ! test the input arguments
           lquery = lwork < 0
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           if (left) then
             lw = n*mb
           else
             lw = m*mb
           end if
           info = 0
           if (.not. left .and. .not. right) then
              info = -1
           else if (.not. tran .and. .not. notran) then
              info = -2
           else if (k < 0) then
             info = -5
           else if (m < k) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < mb .or. mb < 1) then
             info = -6
           else if (lda < max(1,k)) then
             info = -9
           else if (ldt < max(1,mb)) then
             info = -11
           else if (ldc < max(1,m)) then
              info = -13
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -15
           end if
           if (info /= 0) then
             call la_xerbla('WLAMSWLQ',-info)
             work(1) = lw
             return
           else if (lquery) then
             work(1) = lw
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((nb <= k) .or. (nb >= max(m,n,k))) then
             call la_wgemlqt(side,trans,m,n,k,mb,a,lda,t,ldt,c,ldc,work,info)

             return
           end if
           if (left .and. tran) then
               ! multiply q to the last block of c
               kk = mod((m - k), (nb - k))
               ctr = (m - k)/(nb - k)
               if (kk > 0) then
                 ii = m - kk + 1
                 call la_wtpmlqt('L','C',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                           1,1),ldc,c(ii,1),ldc,work,info)
               else
                 ii = m + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+nb)
                 ctr = ctr - 1
                 call la_wtpmlqt('L','C',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                           1),ldc,c(i,1),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:nb)
               call la_wgemlqt('L','C',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (left .and. notran) then
               ! multiply q to the first block of c
              kk = mod((m - k), (nb - k))
              ii = m - kk + 1
              ctr = 1
              call la_wgemlqt('L','N',nb,n,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (i:i+nb,1:n)
               call la_wtpmlqt('L','N',nb - k,n,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(i,1),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= m) then
               ! multiply q to the last block of c
               call la_wtpmlqt('L','N',kk,n,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt, &
                         c(1,1),ldc,c(ii,1),ldc,work,info)
              end if
           else if (right .and. notran) then
               ! multiply q to the last block of c
               kk = mod((n - k), (nb - k))
               ctr = (n - k)/(nb - k)
               if (kk > 0) then
                 ii = n - kk + 1
                 call la_wtpmlqt('R','N',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1), &
                           ldt,c(1,1),ldc,c(1,ii),ldc,work,info)
               else
                 ii = n + 1
               end if
               do i = ii - (nb - k),nb + 1,-(nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               ctr = ctr - 1
               call la_wtpmlqt('R','N',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1), &
                         ldt,c(1,1),ldc,c(1,i),ldc,work,info)
               end do
               ! multiply q to the first block of c (1:m,1:mb)
               call la_wgemlqt('R','N',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                         info)
           else if (right .and. tran) then
             ! multiply q to the first block of c
              kk = mod((n - k), (nb - k))
              ii = n - kk + 1
              call la_wgemlqt('R','C',m,nb,k,mb,a(1,1),lda,t,ldt,c(1,1),ldc,work, &
                        info)
              ctr = 1
              do i = nb + 1,ii - nb + k, (nb - k)
               ! multiply q to the current block of c (1:m,i:i+mb)
               call la_wtpmlqt('R','C',m,nb - k,k,0,mb,a(1,i),lda,t(1,ctr*k + 1),ldt,c(1, &
                         1),ldc,c(1,i),ldc,work,info)
               ctr = ctr + 1
              end do
              if (ii <= n) then
             ! multiply q to the last block of c
               call la_wtpmlqt('R','C',m,kk,k,0,mb,a(1,ii),lda,t(1,ctr*k + 1),ldt,c( &
                         1,1),ldc,c(1,ii),ldc,work,info)
              end if
           end if
           work(1) = lw
           return
     end subroutine la_wlamswlq
#endif

     !> CLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a complex M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_claswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb <= 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('CLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_cgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_cgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_ctplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1),ldt, &
                        work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_ctplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1),ldt, &
                       work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_claswlq
     !> ZLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a complexx M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_zlaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb <= 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('ZLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_zgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_zgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_ztplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_ztplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_zlaswlq
#ifdef LA_WITH_XDP
     !> YLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a complexx M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_ylaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb <= 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('YLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_ygelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_ygelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_ytplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_ytplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_ylaswlq
#endif
#ifdef LA_WITH_QP
     !> WLASWLQ: computes a blocked Tall-Skinny LQ factorization of
     !> a complexx M-by-N matrix A for M <= N:
     !> A = ( L 0 ) *  Q,
     !> where:
     !> Q is a n-by-N orthogonal matrix, stored on exit in an implicit
     !> form in the elements above the diagonal of the array A and in
     !> the elements of the array T;
     !> L is a lower-triangular M-by-M matrix stored on exit in
     !> the elements on and below the diagonal of the array A.
     !> 0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.

     pure subroutine la_wlaswlq(m,n,mb,nb,a,lda,t,ldt,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,mb,nb,lwork,ldt
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: work(*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,ii,kk,ctr
           ! External Subroutines
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
             info = -1
           else if (n < 0 .or. n < m) then
             info = -2
           else if (mb < 1 .or. (mb > m .and. m > 0)) then
             info = -3
           else if (nb <= 0) then
             info = -4
           else if (lda < max(1,m)) then
             info = -6
           else if (ldt < mb) then
             info = -8
           else if ((lwork < m*mb) .and. (.not. lquery)) then
             info = -10
           end if
           if (info == 0) then
           work(1) = mb*m
           end if
           if (info /= 0) then
             call la_xerbla('WLASWLQ',-info)
             return
           else if (lquery) then
            return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
               return
           end if
           ! the lq decomposition
            if ((m >= n) .or. (nb <= m) .or. (nb >= n)) then
             call la_wgelqt(m,n,mb,a,lda,t,ldt,work,info)
             return
            end if
            kk = mod((n - m), (nb - m))
            ii = n - kk + 1
            ! compute the lq factorization of the first block a(1:m,1:nb)
            call la_wgelqt(m,nb,mb,a(1,1),lda,t,ldt,work,info)
            ctr = 1
            do i = nb + 1,ii - nb + m, (nb - m)
            ! compute the qr factorization of the current block a(1:m,i:i+nb-m)
              call la_wtplqt(m,nb - m,0,mb,a(1,1),lda,a(1,i),lda,t(1,ctr*m + 1), &
                        ldt,work,info)
              ctr = ctr + 1
            end do
           ! compute the qr factorization of the last block a(1:m,ii:n)
            if (ii <= n) then
             call la_wtplqt(m,kk,0,mb,a(1,1),lda,a(1,ii),lda,t(1,ctr*m + 1), &
                       ldt,work,info)
            end if
           work(1) = m*mb
           return
     end subroutine la_wlaswlq
#endif

     !> CGELQ: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_cgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'CGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'CGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('CGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_cgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_claswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_cgelq
     !> ZGELQ: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_zgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'ZGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'ZGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('ZGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_zgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_zlaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_zgelq
#ifdef LA_WITH_XDP
     !> YGELQ: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_ygelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'YGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'YGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('YGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_ygelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_ylaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_ygelq
#endif
#ifdef LA_WITH_QP
     !> WGELQ: computes an LQ factorization of a complex M-by-N matrix A:
     !> A = ( L 0 ) *  Q
     !> where:
     !> Q is a N-by-N orthogonal matrix;
     !> L is a lower-triangular M-by-M matrix;
     !> 0 is a M-by-(N-M) zero matrix, if M < N.

     pure subroutine la_wgelq(m,n,a,lda,t,tsize,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd. --
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,tsize,lwork
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: t(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,lminws,mint,minw
           integer(ilp) :: mb,nb,mintsz,nblcks,lwmin,lwopt,lwreq
           ! Intrinsic Functions
           intrinsic :: max,min,mod
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (tsize == -1 .or. tsize == -2 .or. lwork == -1 .or. lwork == -2)
           mint = .false.
           minw = .false.
           if (tsize == -2 .or. lwork == -2) then
             if (tsize /= -1) mint = .true.
             if (lwork /= -1) minw = .true.
           end if
           ! determine the block size
           if (min(m,n) > 0) then
             mb = la_ilaenv(1,'WGELQ ',' ',m,n,1,-1)
             nb = la_ilaenv(1,'WGELQ ',' ',m,n,2,-1)
           else
             mb = 1
             nb = n
           end if
           if (mb > min(m,n) .or. mb < 1) mb = 1
           if (nb > n .or. nb <= m) nb = n
           mintsz = m + 5
           if (nb > m .and. n > m) then
             if (mod(n - m,nb - m) == 0) then
               nblcks = (n - m)/(nb - m)
             else
               nblcks = (n - m)/(nb - m) + 1
             end if
           else
             nblcks = 1
           end if
           ! determine if the workspace size satisfies minimal size
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwmin = max(1,n)
              lwopt = max(1,mb*n)
           else
              lwmin = max(1,m)
              lwopt = max(1,mb*m)
           end if
           lminws = .false.
           if ((tsize < max(1,mb*m*nblcks + 5) .or. lwork < lwopt) .and. (lwork >= lwmin) .and. ( &
                     tsize >= mintsz) .and. (.not. lquery)) then
             if (tsize < max(1,mb*m*nblcks + 5)) then
                 lminws = .true.
                 mb = 1
                 nb = n
             end if
             if (lwork < lwopt) then
                 lminws = .true.
                 mb = 1
             end if
           end if
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
              lwreq = max(1,mb*n)
           else
              lwreq = max(1,mb*m)
           end if
           if (m < 0) then
             info = -1
           else if (n < 0) then
             info = -2
           else if (lda < max(1,m)) then
             info = -4
           else if (tsize < max(1,mb*m*nblcks + 5) .and. (.not. lquery) .and. (.not. lminws)) &
                     then
             info = -6
           else if ((lwork < lwreq) .and. (.not. lquery) .and. (.not. lminws)) then
             info = -8
           end if
           if (info == 0) then
             if (mint) then
               t(1) = mintsz
             else
               t(1) = mb*m*nblcks + 5
             end if
             t(2) = mb
             t(3) = nb
             if (minw) then
               work(1) = lwmin
             else
               work(1) = lwreq
             end if
           end if
           if (info /= 0) then
             call la_xerbla('WGELQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n) == 0) then
             return
           end if
           ! the lq decomposition
           if ((n <= m) .or. (nb <= m) .or. (nb >= n)) then
             call la_wgelqt(m,n,mb,a,lda,t(6),mb,work,info)
           else
             call la_wlaswlq(m,n,mb,nb,a,lda,t(6),mb,work,lwork,info)
           end if
           work(1) = lwreq
           return
     end subroutine la_wgelq
#endif

     !> CGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product
     !> of blocked elementary reflectors computed by short wide
     !> LQ factorization (CGELQ)

     pure subroutine la_cgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),t(*)
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = real(lw,KIND=sp)
           end if
           if (info /= 0) then
             call la_xerbla('CGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_cgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_clamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = real(lw,KIND=sp)
           return
     end subroutine la_cgemlq
     !> ZGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product
     !> of blocked elementary reflectors computed by short wide
     !> LQ factorization (ZGELQ)

     pure subroutine la_zgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),t(*)
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('ZGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_zgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_zlamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_zgemlq
#ifdef LA_WITH_XDP
     !> YGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product
     !> of blocked elementary reflectors computed by short wide
     !> LQ factorization (YGELQ)

     pure subroutine la_ygemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),t(*)
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('YGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_ygemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_ylamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_ygemlq
#endif
#ifdef LA_WITH_QP
     !> WGEMLQ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product
     !> of blocked elementary reflectors computed by short wide
     !> LQ factorization (WGELQ)

     pure subroutine la_wgemlq(side,trans,m,n,k,a,lda,t,tsize,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n,k,tsize,lwork,ldc
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),t(*)
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(out) :: work(*)
       ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,right,tran,notran,lquery
           integer(ilp) :: mb,nb,lw,nblcks,mn
           ! Intrinsic Functions
           intrinsic :: int,max,min,mod
           ! Executable Statements
           ! test the input arguments
           lquery = lwork == -1
           notran = la_lsame(trans,'N')
           tran = la_lsame(trans,'C')
           left = la_lsame(side,'L')
           right = la_lsame(side,'R')
           mb = int(t(2),KIND=ilp)
           nb = int(t(3),KIND=ilp)
           if (left) then
             lw = n*mb
             mn = m
           else
             lw = m*mb
             mn = n
           end if
           if ((nb > k) .and. (mn > k)) then
             if (mod(mn - k,nb - k) == 0) then
               nblcks = (mn - k)/(nb - k)
             else
               nblcks = (mn - k)/(nb - k) + 1
             end if
           else
             nblcks = 1
           end if
           info = 0
           if (.not. left .and. .not. right) then
             info = -1
           else if (.not. tran .and. .not. notran) then
             info = -2
           else if (m < 0) then
             info = -3
           else if (n < 0) then
             info = -4
           else if (k < 0 .or. k > mn) then
             info = -5
           else if (lda < max(1,k)) then
             info = -7
           else if (tsize < 5) then
             info = -9
           else if (ldc < max(1,m)) then
             info = -11
           else if ((lwork < max(1,lw)) .and. (.not. lquery)) then
             info = -13
           end if
           if (info == 0) then
             work(1) = lw
           end if
           if (info /= 0) then
             call la_xerbla('WGEMLQ',-info)
             return
           else if (lquery) then
             return
           end if
           ! quick return if possible
           if (min(m,n,k) == 0) then
             return
           end if
           if ((left .and. m <= k) .or. (right .and. n <= k) .or. (nb <= k) .or. (nb >= max(m,n, &
                     k))) then
             call la_wgemlqt(side,trans,m,n,k,mb,a,lda,t(6),mb,c,ldc,work,info &
                       )
           else
             call la_wlamswlq(side,trans,m,n,k,mb,nb,a,lda,t(6),mb,c,ldc,work, &
                       lwork,info)
           end if
           work(1) = lw
           return
     end subroutine la_wgemlq
#endif

end module la_lapack_orthogonal_factors_ql
