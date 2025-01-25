module la_linalg_lapack_aux
     use la_linalg_constants
     use la_linalg_blas
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_chla_transtype
     public :: la_droundup_lwork
     public :: la_icmax1
     public :: la_ieeeck
     public :: la_ilaclc
     public :: la_ilaclr
     public :: la_iladiag
     public :: la_iladlc
     public :: la_iladlr
     public :: la_ilaenv
     public :: la_ilaenv2stage
     public :: la_ilaprec
     public :: la_ilaslc
     public :: la_ilaslr
     public :: la_ilatrans
     public :: la_ilauplo
     public :: la_ilazlc
     public :: la_ilazlr
     public :: la_iparam2stage
     public :: la_iparmq
     public :: la_izmax1
     public :: la_lsamen
     public :: la_sroundup_lwork
     public :: la_xerbla
     public :: la_xerbla_array
     public :: la_qroundup_lwork
     public :: la_ilaqiag
     public :: la_ilaqlc
     public :: la_ilaqlr
     public :: la_ilawlc
     public :: la_ilawlr
     public :: la_iwmax1
     public :: la_selctg_s
     public :: la_select_s
     public :: la_selctg_d
     public :: la_select_d
     public :: la_selctg_q
     public :: la_select_q
     public :: la_selctg_c
     public :: la_select_c
     public :: la_selctg_z
     public :: la_select_z
     public :: la_selctg_w
     public :: la_select_w

     ! SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments
     ! used to select eigenvalues to sort to the top left of the Schur form.
     ! An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if SELCTG is true, i.e.,
     abstract interface
        pure logical(lk) function la_selctg_s(alphar,alphai,beta)
            import sp,dp,qp,lk
            implicit none
            real(sp),intent(in) :: alphar,alphai,beta
        end function la_selctg_s
        pure logical(lk) function la_select_s(alphar,alphai)
            import sp,dp,qp,lk
            implicit none
            real(sp),intent(in) :: alphar,alphai
        end function la_select_s
        pure logical(lk) function la_selctg_d(alphar,alphai,beta)
            import sp,dp,qp,lk
            implicit none
            real(dp),intent(in) :: alphar,alphai,beta
        end function la_selctg_d
        pure logical(lk) function la_select_d(alphar,alphai)
            import sp,dp,qp,lk
            implicit none
            real(dp),intent(in) :: alphar,alphai
        end function la_select_d
        pure logical(lk) function la_selctg_q(alphar,alphai,beta)
            import sp,dp,qp,lk
            implicit none
            real(qp),intent(in) :: alphar,alphai,beta
        end function la_selctg_q
        pure logical(lk) function la_select_q(alphar,alphai)
            import sp,dp,qp,lk
            implicit none
            real(qp),intent(in) :: alphar,alphai
        end function la_select_q
        pure logical(lk) function la_selctg_c(alpha,beta)
            import sp,dp,qp,lk
            implicit none
            complex(sp),intent(in) :: alpha,beta
        end function la_selctg_c
        pure logical(lk) function la_select_c(alpha)
            import sp,dp,qp,lk
            implicit none
            complex(sp),intent(in) :: alpha
        end function la_select_c
        pure logical(lk) function la_selctg_z(alpha,beta)
            import sp,dp,qp,lk
            implicit none
            complex(dp),intent(in) :: alpha,beta
        end function la_selctg_z
        pure logical(lk) function la_select_z(alpha)
            import sp,dp,qp,lk
            implicit none
            complex(dp),intent(in) :: alpha
        end function la_select_z
        pure logical(lk) function la_selctg_w(alpha,beta)
            import sp,dp,qp,lk
            implicit none
            complex(qp),intent(in) :: alpha,beta
        end function la_selctg_w
        pure logical(lk) function la_select_w(alpha)
            import sp,dp,qp,lk
            implicit none
            complex(qp),intent(in) :: alpha
        end function la_select_w
     end interface

     contains

     !> This subroutine translates from a BLAST-specified integer constant to
     !> the character string specifying a transposition operation.
     !> CHLA_TRANSTYPE: returns an CHARACTER*1.  If CHLA_TRANSTYPE: is 'X',
     !> then input is not an integer indicating a transposition operator.
     !> Otherwise CHLA_TRANSTYPE returns the constant value corresponding to
     !> TRANS.

     pure character function la_chla_transtype(trans)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: trans
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_no_trans = 111
           integer(ilp),parameter :: blas_trans = 112
           integer(ilp),parameter :: blas_conj_trans = 113
           
           ! Executable Statements
           if (trans == blas_no_trans) then
              la_chla_transtype = 'N'
           else if (trans == blas_trans) then
              la_chla_transtype = 'T'
           else if (trans == blas_conj_trans) then
              la_chla_transtype = 'C'
           else
              la_chla_transtype = 'X'
           end if
           return
     end function la_chla_transtype

     !> DROUNDUP_LWORK: deals with a subtle bug with returning LWORK as a Float.
     !> This routine guarantees it is rounded up instead of down by
     !> multiplying LWORK by 1+eps when it is necessary, where eps is the relative machine precision.
     !> E.g.,
     !> float( 9007199254740993            ) == 9007199254740992
     !> float( 9007199254740993 ) * (1.+eps) == 9007199254740994
     !> \return DROUNDUP_LWORK
     !>
     !> DROUNDUP_LWORK >= LWORK.
     !> DROUNDUP_LWORK is guaranteed to have zero decimal part.

     pure real(dp) function la_droundup_lwork(lwork)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lwork
       ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: epsilon,real,int
           ! Executable Statements
           la_droundup_lwork = real(lwork,KIND=dp)
           if (int(la_droundup_lwork,KIND=ilp) < lwork) then
               ! force round up of lwork
               la_droundup_lwork = la_droundup_lwork*(1.0e+0_dp + epsilon(0.0e+0_dp))
                         
           end if
           return
     end function la_droundup_lwork

     !> ICMAX1: finds the index of the first vector element of maximum absolute value.
     !> Based on ICAMAX from Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure integer(ilp) function la_icmax1(n,cx,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: smax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_icmax1 = 0
           if (n < 1 .or. incx <= 0) return
           la_icmax1 = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              smax = abs(cx(1))
              do i = 2,n
                 if (abs(cx(i)) > smax) then
                    la_icmax1 = i
                    smax = abs(cx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              smax = abs(cx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(cx(ix)) > smax) then
                    la_icmax1 = i
                    smax = abs(cx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_icmax1

     !> IEEECK: is called from the ILAENV to verify that Infinity and
     !> possibly NaN arithmetic is safe (i.e. will not trap).

     pure integer(ilp) function la_ieeeck(ispec,zero,one)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ispec
           real(sp),intent(in) :: one,zero
        ! =====================================================================
           ! Local Scalars
           real(sp) :: nan1,nan2,nan3,nan4,nan5,nan6,neginf,negzro,newzro,posinf
           ! Executable Statements
           la_ieeeck = 1
           posinf = one/zero
           if (posinf <= one) then
              la_ieeeck = 0
              return
           end if
           neginf = -one/zero
           if (neginf >= zero) then
              la_ieeeck = 0
              return
           end if
           negzro = one/(neginf + one)
           if (negzro /= zero) then
              la_ieeeck = 0
              return
           end if
           neginf = one/negzro
           if (neginf >= zero) then
              la_ieeeck = 0
              return
           end if
           newzro = negzro + zero
           if (newzro /= zero) then
              la_ieeeck = 0
              return
           end if
           posinf = one/newzro
           if (posinf <= one) then
              la_ieeeck = 0
              return
           end if
           neginf = neginf*posinf
           if (neginf >= zero) then
              la_ieeeck = 0
              return
           end if
           posinf = posinf*posinf
           if (posinf <= one) then
              la_ieeeck = 0
              return
           end if
           ! return if we were only asked to check infinity arithmetic
           if (ispec == 0) return
           nan1 = posinf + neginf
           nan2 = posinf/neginf
           nan3 = posinf/posinf
           nan4 = posinf*zero
           nan5 = neginf*negzro
           nan6 = nan5*zero
           if (nan1 == nan1) then
              la_ieeeck = 0
              return
           end if
           if (nan2 == nan2) then
              la_ieeeck = 0
              return
           end if
           if (nan3 == nan3) then
              la_ieeeck = 0
              return
           end if
           if (nan4 == nan4) then
              la_ieeeck = 0
              return
           end if
           if (nan5 == nan5) then
              la_ieeeck = 0
              return
           end if
           if (nan6 == nan6) then
              la_ieeeck = 0
              return
           end if
           return
     end function la_ieeeck

     !> ILACLC: scans A for its last non-zero column.

     pure integer(ilp) function la_ilaclc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(sp),parameter :: zero = (0.0e+0,0.0e+0)
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_ilaclc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_ilaclc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_ilaclc = n,1,-1
                 do i = 1,m
                    if (a(i,la_ilaclc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_ilaclc

     !> ILACLR: scans A for its last non-zero row.

     pure integer(ilp) function la_ilaclr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(sp),parameter :: zero = (0.0e+0,0.0e+0)
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_ilaclr = m
           else if (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_ilaclr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_ilaclr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_ilaclr = max(la_ilaclr,i)
              end do
           end if
           return
     end function la_ilaclr

     !> This subroutine translated from a character string specifying if a
     !> matrix has unit diagonal or not to the relevant BLAST-specified
     !> integer constant.
     !> ILADIAG: returns an INTEGER.  If ILADIAG: < 0, then the input is not a
     !> character indicating a unit or non-unit diagonal.  Otherwise ILADIAG
     !> returns the constant value corresponding to DIAG.

     integer(ilp) function la_iladiag(diag)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character :: diag
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_non_unit_diag = 131
           integer(ilp),parameter :: blas_unit_diag = 132
           
           ! Executable Statements
           if (la_lsame(diag,'N')) then
              la_iladiag = blas_non_unit_diag
           else if (la_lsame(diag,'U')) then
              la_iladiag = blas_unit_diag
           else
              la_iladiag = -1
           end if
           return
     end function la_iladiag

     !> ILADLC: scans A for its last non-zero column.

     pure integer(ilp) function la_iladlc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: zero = 0.0d+0
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_iladlc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_iladlc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_iladlc = n,1,-1
                 do i = 1,m
                    if (a(i,la_iladlc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_iladlc

     !> ILADLR: scans A for its last non-zero row.

     pure integer(ilp) function la_iladlr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: zero = 0.0d+0
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_iladlr = m
           else if (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_iladlr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_iladlr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_iladlr = max(la_iladlr,i)
              end do
           end if
           return
     end function la_iladlr

     !> This subroutine translated from a character string specifying an
     !> intermediate precision to the relevant BLAST-specified integer
     !> constant.
     !> ILAPREC: returns an INTEGER.  If ILAPREC: < 0, then the input is not a
     !> character indicating a supported intermediate precision.  Otherwise
     !> ILAPREC returns the constant value corresponding to PREC.

     integer(ilp) function la_ilaprec(prec)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character :: prec
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_prec_single = 211
           integer(ilp),parameter :: blas_prec_double = 212
           integer(ilp),parameter :: blas_prec_indigenous = 213
           integer(ilp),parameter :: blas_prec_extra = 214
           
           ! Executable Statements
           if (la_lsame(prec,'S')) then
              la_ilaprec = blas_prec_single
           else if (la_lsame(prec,'D')) then
              la_ilaprec = blas_prec_double
           else if (la_lsame(prec,'I')) then
              la_ilaprec = blas_prec_indigenous
           else if (la_lsame(prec,'X') .or. la_lsame(prec,'E')) then
              la_ilaprec = blas_prec_extra
           else
              la_ilaprec = -1
           end if
           return
     end function la_ilaprec

     !> ILASLC: scans A for its last non-zero column.

     pure integer(ilp) function la_ilaslc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: zero = 0.0e+0
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_ilaslc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_ilaslc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_ilaslc = n,1,-1
                 do i = 1,m
                    if (a(i,la_ilaslc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_ilaslc

     !> ILASLR: scans A for its last non-zero row.

     pure integer(ilp) function la_ilaslr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: zero = 0.0e+0
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_ilaslr = m
           elseif (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_ilaslr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_ilaslr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_ilaslr = max(la_ilaslr,i)
              end do
           end if
           return
     end function la_ilaslr

     !> This subroutine translates from a character string specifying a
     !> transposition operation to the relevant BLAST-specified integer
     !> constant.
     !> ILATRANS: returns an INTEGER.  If ILATRANS: < 0, then the input is not
     !> a character indicating a transposition operator.  Otherwise ILATRANS
     !> returns the constant value corresponding to TRANS.

     integer(ilp) function la_ilatrans(trans)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character :: trans
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_no_trans = 111
           integer(ilp),parameter :: blas_trans = 112
           integer(ilp),parameter :: blas_conj_trans = 113
           
           ! Executable Statements
           if (la_lsame(trans,'N')) then
              la_ilatrans = blas_no_trans
           else if (la_lsame(trans,'T')) then
              la_ilatrans = blas_trans
           else if (la_lsame(trans,'C')) then
              la_ilatrans = blas_conj_trans
           else
              la_ilatrans = -1
           end if
           return
     end function la_ilatrans

     !> This subroutine translated from a character string specifying a
     !> upper- or lower-triangular matrix to the relevant BLAST-specified
     !> integer constant.
     !> ILAUPLO: returns an INTEGER.  If ILAUPLO: < 0, then the input is not
     !> a character indicating an upper- or lower-triangular matrix.
     !> Otherwise ILAUPLO returns the constant value corresponding to UPLO.

     integer(ilp) function la_ilauplo(uplo)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character :: uplo
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_upper = 121
           integer(ilp),parameter :: blas_lower = 122
           
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              la_ilauplo = blas_upper
           else if (la_lsame(uplo,'L')) then
              la_ilauplo = blas_lower
           else
              la_ilauplo = -1
           end if
           return
     end function la_ilauplo

     !> ILAZLC: scans A for its last non-zero column.

     pure integer(ilp) function la_ilazlc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(dp),parameter :: zero = (0.0d+0,0.0d+0)
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_ilazlc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_ilazlc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_ilazlc = n,1,-1
                 do i = 1,m
                    if (a(i,la_ilazlc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_ilazlc

     !> ILAZLR: scans A for its last non-zero row.

     pure integer(ilp) function la_ilazlr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(dp),parameter :: zero = (0.0d+0,0.0d+0)
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_ilazlr = m
           else if (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_ilazlr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_ilazlr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_ilazlr = max(la_ilazlr,i)
              end do
           end if
           return
     end function la_ilazlr

     !> This program sets problem and machine dependent parameters
     !> useful for xHSEQR and related subroutines for eigenvalue
     !> problems. It is called whenever
     !> IPARMQ: is called with 12 <= ISPEC <= 16

     pure integer(ilp) function la_iparmq(ispec,name,opts,n,ilo,ihi,lwork)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ihi,ilo,ispec,lwork,n
           character,intent(in) :: name*(*),opts*(*)
        ! ================================================================
           ! Parameters
           integer(ilp),parameter :: inmin = 12
           integer(ilp),parameter :: inwin = 13
           integer(ilp),parameter :: inibl = 14
           integer(ilp),parameter :: ishfts = 15
           integer(ilp),parameter :: iacc22 = 16
           integer(ilp),parameter :: icost = 17
           integer(ilp),parameter :: nmin = 75
           integer(ilp),parameter :: k22min = 14
           integer(ilp),parameter :: kacmin = 14
           integer(ilp),parameter :: nibble = 14
           integer(ilp),parameter :: knwswp = 500
           integer(ilp),parameter :: rcost = 10
           real(sp),parameter :: two = 2.0
           
           ! Local Scalars
           integer(ilp) :: nh,ns
           integer(ilp) :: i,ic,iz
           character :: subnam*6
           ! Intrinsic Functions
           intrinsic :: log,max,mod,nint,real
           ! Executable Statements
           if ((ispec == ishfts) .or. (ispec == inwin) .or. (ispec == iacc22)) then
              ! ==== set the number simultaneous shifts ====
              nh = ihi - ilo + 1
              ns = 2
              if (nh >= 30) ns = 4
              if (nh >= 60) ns = 10
              if (nh >= 150) ns = max(10,nh/nint(log(real(nh,KIND=dp))/log(two), &
                        KIND=ilp))
              if (nh >= 590) ns = 64
              if (nh >= 3000) ns = 128
              if (nh >= 6000) ns = 256
              ns = max(2,ns - mod(ns,2))
           end if
           if (ispec == inmin) then
              ! ===== matrices of order smaller than nmin get sent
              ! .     to xlahqr, the classic double shift algorithm.
              ! .     this must be at least 11. ====
              la_iparmq = nmin
           else if (ispec == inibl) then
              ! ==== inibl: skip a multi-shift qr iteration and
              ! .    whenever aggressive early deflation finds
              ! .    at least (nibble*(window size)/100) deflations. ====
              la_iparmq = nibble
           else if (ispec == ishfts) then
              ! ==== nshfts: the number of simultaneous shifts =====
              la_iparmq = ns
           else if (ispec == inwin) then
              ! ==== nw: deflation window size.  ====
              if (nh <= knwswp) then
                 la_iparmq = ns
              else
                 la_iparmq = 3*ns/2
              end if
           else if (ispec == iacc22) then
              ! ==== iacc22: whether to accumulate reflections
              ! .     before updating the far-from-diagonal elements
              ! .     and whether to use 2-by-2 block structure while
              ! .     doing it.  a small amount of work could be saved
              ! .     by making this choice dependent also upon the
              ! .     nh=ihi-ilo+1.
              ! convert name to upper case if the first character is lower case.
              la_iparmq = 0
              subnam = name
              ic = ichar(subnam(1:1))
              iz = ichar('Z')
              if (iz == 90 .or. iz == 122) then
                 ! ascii character set
                 if (ic >= 97 .and. ic <= 122) then
                    subnam(1:1) = char(ic - 32)
                    do i = 2,6
                       ic = ichar(subnam(i:i))
                       if (ic >= 97 .and. ic <= 122) subnam(i:i) = char(ic - 32)
                    end do
                 end if
              else if (iz == 233 .or. iz == 169) then
                 ! ebcdic character set
                 if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 .and. &
                           ic <= 169)) then
                    subnam(1:1) = char(ic + 64)
                    do i = 2,6
                       ic = ichar(subnam(i:i))
                       if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 &
                                 .and. ic <= 169)) subnam(i:i) = char(ic + 64)
                    end do
                 end if
              else if (iz == 218 .or. iz == 250) then
                 ! prime machines:  ascii+128
                 if (ic >= 225 .and. ic <= 250) then
                    subnam(1:1) = char(ic - 32)
                    do i = 2,6
                       ic = ichar(subnam(i:i))
                       if (ic >= 225 .and. ic <= 250) subnam(i:i) = char(ic - 32)
                    end do
                 end if
              end if
              if (subnam(2:6) == 'GGHRD' .or. subnam(2:6) == 'GGHD3') then
                 la_iparmq = 1
                 if (nh >= k22min) la_iparmq = 2
              else if (subnam(4:6) == 'EXC') then
                 if (nh >= kacmin) la_iparmq = 1
                 if (nh >= k22min) la_iparmq = 2
              else if (subnam(2:6) == 'HSEQR' .or. subnam(2:5) == 'LAQR') then
                 if (ns >= kacmin) la_iparmq = 1
                 if (ns >= k22min) la_iparmq = 2
              end if
           else if (ispec == icost) then
              ! === relative cost of near-the-diagonal chase vs
                  ! blas updates ===
              la_iparmq = rcost
           else
              ! ===== invalid value of ispec =====
              la_iparmq = -1
           end if
     end function la_iparmq

     !> IZMAX1: finds the index of the first vector element of maximum absolute value.
     !> Based on IZAMAX from Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure integer(ilp) function la_izmax1(n,zx,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dmax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_izmax1 = 0
           if (n < 1 .or. incx <= 0) return
           la_izmax1 = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = abs(zx(1))
              do i = 2,n
                 if (abs(zx(i)) > dmax) then
                    la_izmax1 = i
                    dmax = abs(zx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = abs(zx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(zx(ix)) > dmax) then
                    la_izmax1 = i
                    dmax = abs(zx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_izmax1

     !> LSAMEN:  tests if the first N letters of CA are the same as the
     !> first N letters of CB, regardless of case.
     !> LSAMEN returns .TRUE. if CA and CB are equivalent except for case
     !> and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
     !> or LEN( CB ) is less than N.

     pure logical(lk) function la_lsamen(n,ca,cb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character(len=*),intent(in) :: ca,cb
           integer(ilp),intent(in) :: n
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i
           ! Intrinsic Functions
           intrinsic :: len
           ! Executable Statements
           la_lsamen = .false.
           if (len(ca) < n .or. len(cb) < n) go to 20
           ! do for each character in the two strings.
           do i = 1,n
              ! test if the characters are equal using la_lsame.
              if (.not. la_lsame(ca(i:i),cb(i:i))) go to 20
           end do
           la_lsamen = .true.
20         continue
           return
     end function la_lsamen

     !> SROUNDUP_LWORK: deals with a subtle bug with returning LWORK as a Float.
     !> This routine guarantees it is rounded up instead of down by
     !> multiplying LWORK by 1+eps when it is necessary, where eps is the relative machine precision.
     !> E.g.,
     !> float( 16777217            ) == 16777216
     !> float( 16777217 ) * (1.+eps) == 16777218
     !> \return SROUNDUP_LWORK
     !>
     !> SROUNDUP_LWORK >= LWORK.
     !> SROUNDUP_LWORK is guaranteed to have zero decimal part.

     pure real(sp) function la_sroundup_lwork(lwork)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lwork
       ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: epsilon,real,int
           ! Executable Statements
           la_sroundup_lwork = real(lwork,KIND=sp)
           if (int(la_sroundup_lwork,KIND=ilp) < lwork) then
               ! force round up of lwork
               la_sroundup_lwork = la_sroundup_lwork*(1.0e+0_sp + epsilon(0.0e+0_sp))
                         
           end if
           return
     end function la_sroundup_lwork

     !> DROUNDUP_LWORK: deals with a subtle bug with returning LWORK as a Float.
     !> This routine guarantees it is rounded up instead of down by
     !> multiplying LWORK by 1+eps when it is necessary, where eps is the relative machine precision.
     !> E.g.,
     !> float( 9007199254740993            ) == 9007199254740992
     !> float( 9007199254740993 ) * (1.+eps) == 9007199254740994
     !> \return DROUNDUP_LWORK
     !>
     !> DROUNDUP_LWORK >= LWORK.
     !> DROUNDUP_LWORK is guaranteed to have zero decimal part.

     pure real(qp) function la_qroundup_lwork(lwork)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lwork
       ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: epsilon,real,int
           ! Executable Statements
           la_qroundup_lwork = real(lwork,KIND=qp)
           if (int(la_qroundup_lwork,KIND=ilp) < lwork) then
               ! force round up of lwork
               la_qroundup_lwork = la_qroundup_lwork*(1.0e+0_qp + epsilon(0.0e+0_qp))
                         
           end if
           return
     end function la_qroundup_lwork

     !> This subroutine translated from a character string specifying if a
     !> matrix has unit diagonal or not to the relevant BLAST-specified
     !> integer constant.
     !> ILADIAG: returns an INTEGER.  If ILADIAG: < 0, then the input is not a
     !> character indicating a unit or non-unit diagonal.  Otherwise ILADIAG
     !> returns the constant value corresponding to DIAG.

     integer(ilp) function la_ilaqiag(diag)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character :: diag
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blas_non_unit_qiag = 131
           integer(ilp),parameter :: blas_unit_qiag = 132
           
           ! Executable Statements
           if (la_lsame(diag,'N')) then
              la_ilaqiag = blas_non_unit_qiag
           else if (la_lsame(diag,'U')) then
              la_ilaqiag = blas_unit_qiag
           else
              la_ilaqiag = -1
           end if
           return
     end function la_ilaqiag

     !> ILADLC: scans A for its last non-zero column.

     pure integer(ilp) function la_ilaqlc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: zero = 0.0d+0
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_ilaqlc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_ilaqlc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_ilaqlc = n,1,-1
                 do i = 1,m
                    if (a(i,la_ilaqlc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_ilaqlc

     !> ILADLR: scans A for its last non-zero row.

     pure integer(ilp) function la_ilaqlr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: zero = 0.0d+0
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_ilaqlr = m
           else if (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_ilaqlr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_ilaqlr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_ilaqlr = max(la_ilaqlr,i)
              end do
           end if
           return
     end function la_ilaqlr

     !> ILAZLC: scans A for its last non-zero column.

     pure integer(ilp) function la_ilawlc(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(qp),parameter :: zero = (0.0d+0,0.0d+0)
           
           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (n == 0) then
              la_ilawlc = n
           else if (a(1,n) /= zero .or. a(m,n) /= zero) then
              la_ilawlc = n
           else
           ! now scan each column from the end, returning with the first non-zero.
              do la_ilawlc = n,1,-1
                 do i = 1,m
                    if (a(i,la_ilawlc) /= zero) return
                 end do
              end do
           end if
           return
     end function la_ilawlc

     !> ILAZLR: scans A for its last non-zero row.

     pure integer(ilp) function la_ilawlr(m,n,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: m,n,lda
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           complex(qp),parameter :: zero = (0.0d+0,0.0d+0)
           
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick test for the common case where one corner is non-zero.
           if (m == 0) then
              la_ilawlr = m
           else if (a(m,1) /= zero .or. a(m,n) /= zero) then
              la_ilawlr = m
           else
           ! scan up each column tracking the last zero row seen.
              la_ilawlr = 0
              do j = 1,n
                 i = m
                 do while ((a(max(i,1),j) == zero) .and. (i >= 1))
                    i = i - 1
                 end do
                 la_ilawlr = max(la_ilawlr,i)
              end do
           end if
           return
     end function la_ilawlr

     !> IZMAX1: finds the index of the first vector element of maximum absolute value.
     !> Based on IZAMAX from Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure integer(ilp) function la_iwmax1(n,zx,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(in) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: dmax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_iwmax1 = 0
           if (n < 1 .or. incx <= 0) return
           la_iwmax1 = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = abs(zx(1))
              do i = 2,n
                 if (abs(zx(i)) > dmax) then
                    la_iwmax1 = i
                    dmax = abs(zx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = abs(zx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(zx(ix)) > dmax) then
                    la_iwmax1 = i
                    dmax = abs(zx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_iwmax1

     !> ILAENV: is called from the LAPACK routines to choose problem-dependent
     !> parameters for the local environment.  See ISPEC for a description of
     !> the parameters.
     !> ILAENV returns an INTEGER
     !> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
     !> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
     !> This version provides a set of parameters which should give good,
     !> but not optimal, performance on many of the currently available
     !> computers.  Users are encouraged to modify this subroutine to set
     !> the tuning parameters for their particular machine using the option
     !> and problem size information in the arguments.
     !> This routine will not function correctly if it is converted to all
     !> lower case.  Converting it to all upper case is allowed.

     pure integer(ilp) function la_ilaenv(ispec,name,opts,n1,n2,n3,n4)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character(len=*),intent(in) :: name,opts
           integer(ilp),intent(in) :: ispec,n1,n2,n3,n4
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,iz,nb,nbmin,nx
           logical(lk) :: cname,sname,twostage
           character :: c1*1,c2*2,c4*2,c3*3,subnam*16
           ! Intrinsic Functions
           intrinsic :: char,ichar,int,min,real
           ! Executable Statements
           go to(10,10,10,80,90,100,110,120,130,140,150,160,160,160,160,160,160) &
                     ispec
           ! invalid value for ispec
           la_ilaenv = -1
           return
10         continue
           ! convert name to upper case if the first character is lower case.
           la_ilaenv = 1
           subnam = name
           ic = ichar(subnam(1:1))
           iz = ichar('Z')
           if (iz == 90 .or. iz == 122) then
              ! ascii character set
              if (ic >= 97 .and. ic <= 122) then
                 subnam(1:1) = char(ic - 32)
                 do i = 2,6
                    ic = ichar(subnam(i:i))
                    if (ic >= 97 .and. ic <= 122) subnam(i:i) = char(ic - 32)
                 end do
              end if
           else if (iz == 233 .or. iz == 169) then
              ! ebcdic character set
              if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 .and. &
                        ic <= 169)) then
                 subnam(1:1) = char(ic + 64)
                 do i = 2,6
                    ic = ichar(subnam(i:i))
                    if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 &
                              .and. ic <= 169)) subnam(i:i) = char(ic + 64)
                 end do
              end if
           else if (iz == 218 .or. iz == 250) then
              ! prime machines:  ascii+128
              if (ic >= 225 .and. ic <= 250) then
                 subnam(1:1) = char(ic - 32)
                 do i = 2,6
                    ic = ichar(subnam(i:i))
                    if (ic >= 225 .and. ic <= 250) subnam(i:i) = char(ic - 32)
                 end do
              end if
           end if
           c1 = subnam(1:1)
           sname = c1 == 'S' .or. c1 == 'D'
           cname = c1 == 'C' .or. c1 == 'Z'
           if (.not. (cname .or. sname)) return
           c2 = subnam(2:3)
           c3 = subnam(4:6)
           c4 = c3(2:3)
           twostage = len(subnam) >= 11 .and. subnam(11:11) == '2'
           go to(50,60,70) ispec
50         continue
           ! ispec = 1:  block size
           ! in these examples, separate code is provided for setting nb for
           ! real and complex.  we assume that nb will take the same value in
           ! single or double precision.
           nb = 1
           if (subnam(2:6) == 'LAORH') then
              ! this is for *laorhr_getrfnp routine
              if (sname) then
                  nb = 32
              else
                  nb = 32
              end if
           else if (c2 == 'GE') then
              if (c3 == 'TRF') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              else if (c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. c3 == 'QLF') then
                 if (sname) then
                    nb = 32
                 else
                    nb = 32
                 end if
              else if (c3 == 'QR ') then
                 if (n3 == 1) then
                    if (sname) then
           ! m*n
                       if ((n1*n2 <= 131072) .or. (n1 <= 8192)) then
                          nb = n1
                       else
                          nb = 32768/n2
                       end if
                    else
                       if ((n1*n2 <= 131072) .or. (n1 <= 8192)) then
                          nb = n1
                       else
                          nb = 32768/n2
                       end if
                    end if
                 else
                    if (sname) then
                       nb = 1
                    else
                       nb = 1
                    end if
                 end if
              else if (c3 == 'LQ ') then
                 if (n3 == 2) then
                    if (sname) then
           ! m*n
                       if ((n1*n2 <= 131072) .or. (n1 <= 8192)) then
                          nb = n1
                       else
                          nb = 32768/n2
                       end if
                    else
                       if ((n1*n2 <= 131072) .or. (n1 <= 8192)) then
                          nb = n1
                       else
                          nb = 32768/n2
                       end if
                    end if
                 else
                    if (sname) then
                       nb = 1
                    else
                       nb = 1
                    end if
                 end if
              else if (c3 == 'HRD') then
                 if (sname) then
                    nb = 32
                 else
                    nb = 32
                 end if
              else if (c3 == 'BRD') then
                 if (sname) then
                    nb = 32
                 else
                    nb = 32
                 end if
              else if (c3 == 'TRI') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              end if
           else if (c2 == 'PO') then
              if (c3 == 'TRF') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              end if
           else if (c2 == 'SY') then
              if (c3 == 'TRF') then
                 if (sname) then
                    if (twostage) then
                       nb = 192
                    else
                       nb = 64
                    end if
                 else
                    if (twostage) then
                       nb = 192
                    else
                       nb = 64
                    end if
                 end if
              else if (sname .and. c3 == 'TRD') then
                 nb = 32
              else if (sname .and. c3 == 'GST') then
                 nb = 64
              end if
           else if (cname .and. c2 == 'HE') then
              if (c3 == 'TRF') then
                 if (twostage) then
                    nb = 192
                 else
                    nb = 64
                 end if
              else if (c3 == 'TRD') then
                 nb = 32
              else if (c3 == 'GST') then
                 nb = 64
              end if
           else if (sname .and. c2 == 'OR') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nb = 32
                 end if
              else if (c3(1:1) == 'M') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nb = 32
                 end if
              end if
           else if (cname .and. c2 == 'UN') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nb = 32
                 end if
              else if (c3(1:1) == 'M') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nb = 32
                 end if
              end if
           else if (c2 == 'GB') then
              if (c3 == 'TRF') then
                 if (sname) then
                    if (n4 <= 64) then
                       nb = 1
                    else
                       nb = 32
                    end if
                 else
                    if (n4 <= 64) then
                       nb = 1
                    else
                       nb = 32
                    end if
                 end if
              end if
           else if (c2 == 'PB') then
              if (c3 == 'TRF') then
                 if (sname) then
                    if (n2 <= 64) then
                       nb = 1
                    else
                       nb = 32
                    end if
                 else
                    if (n2 <= 64) then
                       nb = 1
                    else
                       nb = 32
                    end if
                 end if
              end if
           else if (c2 == 'TR') then
              if (c3 == 'TRI') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              else if (c3 == 'EVC') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              end if
           else if (c2 == 'LA') then
              if (c3 == 'UUM') then
                 if (sname) then
                    nb = 64
                 else
                    nb = 64
                 end if
              end if
           else if (sname .and. c2 == 'ST') then
              if (c3 == 'EBZ') then
                 nb = 1
              end if
           else if (c2 == 'GG') then
              nb = 32
              if (c3 == 'HD3') then
                 if (sname) then
                    nb = 32
                 else
                    nb = 32
                 end if
              end if
           end if
           la_ilaenv = nb
           return
60         continue
           ! ispec = 2:  minimum block size
           nbmin = 2
           if (c2 == 'GE') then
              if (c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. c3 == 'QLF') then
                 if (sname) then
                    nbmin = 2
                 else
                    nbmin = 2
                 end if
              else if (c3 == 'HRD') then
                 if (sname) then
                    nbmin = 2
                 else
                    nbmin = 2
                 end if
              else if (c3 == 'BRD') then
                 if (sname) then
                    nbmin = 2
                 else
                    nbmin = 2
                 end if
              else if (c3 == 'TRI') then
                 if (sname) then
                    nbmin = 2
                 else
                    nbmin = 2
                 end if
              end if
           else if (c2 == 'SY') then
              if (c3 == 'TRF') then
                 if (sname) then
                    nbmin = 8
                 else
                    nbmin = 8
                 end if
              else if (sname .and. c3 == 'TRD') then
                 nbmin = 2
              end if
           else if (cname .and. c2 == 'HE') then
              if (c3 == 'TRD') then
                 nbmin = 2
              end if
           else if (sname .and. c2 == 'OR') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nbmin = 2
                 end if
              else if (c3(1:1) == 'M') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nbmin = 2
                 end if
              end if
           else if (cname .and. c2 == 'UN') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nbmin = 2
                 end if
              else if (c3(1:1) == 'M') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nbmin = 2
                 end if
              end if
           else if (c2 == 'GG') then
              nbmin = 2
              if (c3 == 'HD3') then
                 nbmin = 2
              end if
           end if
           la_ilaenv = nbmin
           return
70         continue
           ! ispec = 3:  crossover point
           nx = 0
           if (c2 == 'GE') then
              if (c3 == 'QRF' .or. c3 == 'RQF' .or. c3 == 'LQF' .or. c3 == 'QLF') then
                 if (sname) then
                    nx = 128
                 else
                    nx = 128
                 end if
              else if (c3 == 'HRD') then
                 if (sname) then
                    nx = 128
                 else
                    nx = 128
                 end if
              else if (c3 == 'BRD') then
                 if (sname) then
                    nx = 128
                 else
                    nx = 128
                 end if
              end if
           else if (c2 == 'SY') then
              if (sname .and. c3 == 'TRD') then
                 nx = 32
              end if
           else if (cname .and. c2 == 'HE') then
              if (c3 == 'TRD') then
                 nx = 32
              end if
           else if (sname .and. c2 == 'OR') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nx = 128
                 end if
              end if
           else if (cname .and. c2 == 'UN') then
              if (c3(1:1) == 'G') then
                 if (c4 == 'QR' .or. c4 == 'RQ' .or. c4 == 'LQ' .or. c4 == 'QL' .or. c4 == 'HR' .or. &
                           c4 == 'TR' .or. c4 == 'BR') then
                    nx = 128
                 end if
              end if
           else if (c2 == 'GG') then
              nx = 128
              if (c3 == 'HD3') then
                 nx = 128
              end if
           end if
           la_ilaenv = nx
           return
80         continue
           ! ispec = 4:  number of shifts (used by xhseqr)
           la_ilaenv = 6
           return
90         continue
           ! ispec = 5:  minimum column dimension (not used)
           la_ilaenv = 2
           return
100        continue
           ! ispec = 6:  crossover point for svd (used by xgelss and xgesvd)
           la_ilaenv = int(real(min(n1,n2),KIND=dp)*1.6e0,KIND=ilp)
           return
110        continue
           ! ispec = 7:  number of processors (not used)
           la_ilaenv = 1
           return
120        continue
           ! ispec = 8:  crossover point for multishift (used by xhseqr)
           la_ilaenv = 50
           return
130        continue
           ! ispec = 9:  maximum size of the subproblems at the bottom of the
                       ! computation tree in the divide-and-conquer algorithm
                       ! (used by xgelsd and xgesdd)
           la_ilaenv = 25
           return
140        continue
           ! ispec = 10: ieee and infinity nan arithmetic can be trusted not to trap
           ! la_ilaenv = 0
           la_ilaenv = 1
           if (la_ilaenv == 1) then
              la_ilaenv = la_ieeeck(1,0.0,1.0)
           end if
           return
150        continue
           ! ispec = 11: ieee infinity arithmetic can be trusted not to trap
           ! la_ilaenv = 0
           la_ilaenv = 1
           if (la_ilaenv == 1) then
              la_ilaenv = la_ieeeck(0,0.0,1.0)
           end if
           return
160        continue
           ! 12 <= ispec <= 17: xhseqr or related subroutines.
           la_ilaenv = la_iparmq(ispec,name,opts,n1,n2,n3,n4)
           return
     end function la_ilaenv

     !> This program sets problem and machine dependent parameters
     !> useful for xHETRD_2STAGE, xHETRD_HE2HB, xHETRD_HB2ST,
     !> xGEBRD_2STAGE, xGEBRD_GE2GB, xGEBRD_GB2BD
     !> and related subroutines for eigenvalue problems.
     !> It is called whenever ILAENV is called with 17 <= ISPEC <= 21.
     !> It is called whenever ILAENV2STAGE is called with 1 <= ISPEC <= 5
     !> with a direct conversion ISPEC + 16.

     pure integer(ilp) function la_iparam2stage(ispec,name,opts,ni,nbi,ibi,nxi)
#if defined(_OPENMP)
#endif
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character(len=*),intent(in) :: name,opts
           integer(ilp),intent(in) :: ispec,ni,nbi,ibi,nxi
        ! ================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,iz,kd,ib,lhous,lwork,nthreads,factoptnb,qroptnb, &
                     lqoptnb
           logical(lk) :: rprec,cprec
           character :: prec*1,algo*3,stag*5,subnam*12,vect*1
           ! Intrinsic Functions
           intrinsic :: char,ichar,max
           ! Executable Statements
           ! invalid value for ispec
           if ((ispec < 17) .or. (ispec > 21)) then
               la_iparam2stage = -1
               return
           end if
           ! get the number of threads
           nthreads = 1
#if defined(_OPENMP)
!$OMP PARALLEL
           nthreads = omp_get_num_threads()
!$OMP END PARALLEL
#endif
            ! write(*,*) 'iparam voici nthreads ispec ',nthreads, ispec
           if (ispec /= 19) then
              ! convert name to upper case if the first character is lower case.
              la_iparam2stage = -1
              subnam = name
              ic = ichar(subnam(1:1))
              iz = ichar('Z')
              if (iz == 90 .or. iz == 122) then
                 ! ascii character set
                 if (ic >= 97 .and. ic <= 122) then
                    subnam(1:1) = char(ic - 32)
                    do i = 2,12
                       ic = ichar(subnam(i:i))
                       if (ic >= 97 .and. ic <= 122) subnam(i:i) = char(ic - 32)
                    end do
                 end if
              else if (iz == 233 .or. iz == 169) then
                 ! ebcdic character set
                 if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 .and. &
                           ic <= 169)) then
                    subnam(1:1) = char(ic + 64)
                    do i = 2,12
                       ic = ichar(subnam(i:i))
                       if ((ic >= 129 .and. ic <= 137) .or. (ic >= 145 .and. ic <= 153) .or. (ic >= 162 &
                                 .and. ic <= 169)) subnam(i:i) = char(ic + 64)
                    end do
                 end if
              else if (iz == 218 .or. iz == 250) then
                 ! prime machines:  ascii+128
                 if (ic >= 225 .and. ic <= 250) then
                    subnam(1:1) = char(ic - 32)
                    do i = 2,12
                      ic = ichar(subnam(i:i))
                      if (ic >= 225 .and. ic <= 250) subnam(i:i) = char(ic - 32)
                    end do
                 end if
              end if
              prec = subnam(1:1)
              algo = subnam(4:6)
              stag = subnam(8:12)
              rprec = prec == 'S' .or. prec == 'D'
              cprec = prec == 'C' .or. prec == 'Z'
              ! invalid value for precision
              if (.not. (rprec .or. cprec)) then
                  la_iparam2stage = -1
                  return
              end if
           end if
            ! write(*,*),'rprec,cprec ',rprec,cprec,
           ! $           '   algo ',algo,'    stage ',stag
           if ((ispec == 17) .or. (ispec == 18)) then
           ! ispec = 17, 18:  block size kd, ib
           ! could be also dependent from n but for now it
           ! depend only on sequential or parallel
              if (nthreads > 4) then
                 if (cprec) then
                    kd = 128
                    ib = 32
                 else
                    kd = 160
                    ib = 40
                 end if
              else if (nthreads > 1) then
                 if (cprec) then
                    kd = 64
                    ib = 32
                 else
                    kd = 64
                    ib = 32
                 end if
              else
                 if (cprec) then
                    kd = 16
                    ib = 16
                 else
                    kd = 32
                    ib = 16
                 end if
              end if
              if (ispec == 17) la_iparam2stage = kd
              if (ispec == 18) la_iparam2stage = ib
           else if (ispec == 19) then
           ! ispec = 19:
           ! lhous length of the houselholder representation
           ! matrix (v,t) of the second stage. should be >= 1.
           ! will add the vect option here next release
              vect = opts(1:1)
              if (vect == 'N') then
                 lhous = max(1,4*ni)
              else
                 ! this is not correct, it need to call the algo and the stage2
                 lhous = max(1,4*ni) + ibi
              end if
              if (lhous >= 0) then
                 la_iparam2stage = lhous
              else
                 la_iparam2stage = -1
              end if
           else if (ispec == 20) then
           ! ispec = 20: (21 for future use)
           ! lwork length of the workspace for
           ! either or both stages for trd and brd. should be >= 1.
           ! trd:
           ! trd_stage 1: = lt + lw + ls1 + ls2
                        ! = ldt*kd + n*kd + n*max(kd,factoptnb) + lds2*kd
                          ! where ldt=lds2=kd
                        ! = n*kd + n*max(kd,factoptnb) + 2*kd*kd
           ! trd_stage 2: = (2nb+1)*n + kd*nthreads
           ! trd_both   : = max(stage1,stage2) + ab ( ab=(kd+1)*n )
                        ! = n*kd + n*max(kd+1,factoptnb)
                          ! + max(2*kd*kd, kd*nthreads)
                          ! + (kd+1)*n
              lwork = -1
              subnam(1:1) = prec
              subnam(2:6) = 'GEQRF'
              qroptnb = la_ilaenv(1,subnam,' ',ni,nbi,-1,-1)
              subnam(2:6) = 'GELQF'
              lqoptnb = la_ilaenv(1,subnam,' ',nbi,ni,-1,-1)
              ! could be qr or lq for trd and the max for brd
              factoptnb = max(qroptnb,lqoptnb)
              if (algo == 'TRD') then
                 if (stag == '2STAG') then
                    lwork = ni*nbi + ni*max(nbi + 1,factoptnb) + max(2*nbi*nbi,nbi*nthreads) + (nbi + &
                              1)*ni
                 else if ((stag == 'HE2HB') .or. (stag == 'SY2SB')) then
                    lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
                 else if ((stag == 'HB2ST') .or. (stag == 'SB2ST')) then
                    lwork = (2*nbi + 1)*ni + nbi*nthreads
                 end if
              else if (algo == 'BRD') then
                 if (stag == '2STAG') then
                    lwork = 2*ni*nbi + ni*max(nbi + 1,factoptnb) + max(2*nbi*nbi,nbi*nthreads) + ( &
                              nbi + 1)*ni
                 else if (stag == 'GE2GB') then
                    lwork = ni*nbi + ni*max(nbi,factoptnb) + 2*nbi*nbi
                 else if (stag == 'GB2BD') then
                    lwork = (3*nbi + 1)*ni + nbi*nthreads
                 end if
              end if
              lwork = max(1,lwork)
              if (lwork > 0) then
                 la_iparam2stage = lwork
              else
                 la_iparam2stage = -1
              end if
           else if (ispec == 21) then
           ! ispec = 21 for future use
              la_iparam2stage = nxi
           end if
     end function la_iparam2stage

     !> ILAENV2STAGE: is called from the LAPACK routines to choose problem-dependent
     !> parameters for the local environment.  See ISPEC for a description of
     !> the parameters.
     !> It sets problem and machine dependent parameters useful for *_2STAGE and
     !> related subroutines.
     !> ILAENV2STAGE returns an INTEGER
     !> if ILAENV2STAGE >= 0: ILAENV2STAGE returns the value of the parameter
     !> specified by ISPEC
     !> if ILAENV2STAGE < 0:  if ILAENV2STAGE = -k, the k-th argument had an
     !> illegal value.
     !> This version provides a set of parameters which should give good,
     !> but not optimal, performance on many of the currently available
     !> computers for the 2-stage solvers. Users are encouraged to modify this
     !> subroutine to set the tuning parameters for their particular machine using
     !> the option and problem size information in the arguments.
     !> This routine will not function correctly if it is converted to all
     !> lower case.  Converting it to all upper case is allowed.

     pure integer(ilp) function la_ilaenv2stage(ispec,name,opts,n1,n2,n3,n4)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! july 2017
           ! Scalar Arguments
           character(len=*),intent(in) :: name,opts
           integer(ilp),intent(in) :: ispec,n1,n2,n3,n4
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: iispec
           ! Executable Statements
           go to(10,10,10,10,10) ispec
           ! invalid value for ispec
           la_ilaenv2stage = -1
           return
10         continue
           ! 2stage eigenvalues and svd or related subroutines.
           iispec = 16 + ispec
           la_ilaenv2stage = la_iparam2stage(iispec,name,opts,n1,n2,n3,n4)
           return
     end function la_ilaenv2stage

end module la_linalg_lapack_aux
