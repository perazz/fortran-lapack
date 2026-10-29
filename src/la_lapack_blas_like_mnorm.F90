!> BLAS-like matrix norms
module la_lapack_blas_like_mnorm
     use la_constants
     use la_blas_aux
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_scalar
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slangb
     public :: la_slangb_sub
     public :: la_slange
     public :: la_slange_sub
     public :: la_slangt
     public :: la_slanhs
     public :: la_slanhs_sub
     public :: la_slansb
     public :: la_slansb_sub
     public :: la_slansf
     public :: la_slansf_sub
     public :: la_slansp
     public :: la_slansp_sub
     public :: la_slanst
     public :: la_slansy
     public :: la_slansy_sub
     public :: la_slantb
     public :: la_slantb_sub
     public :: la_slantp
     public :: la_slantp_sub
     public :: la_slantr
     public :: la_slantr_sub
     public :: la_dlangb
     public :: la_dlangb_sub
     public :: la_dlange
     public :: la_dlange_sub
     public :: la_dlangt
     public :: la_dlanhs
     public :: la_dlanhs_sub
     public :: la_dlansb
     public :: la_dlansb_sub
     public :: la_dlansf
     public :: la_dlansf_sub
     public :: la_dlansp
     public :: la_dlansp_sub
     public :: la_dlanst
     public :: la_dlansy
     public :: la_dlansy_sub
     public :: la_dlantb
     public :: la_dlantb_sub
     public :: la_dlantp
     public :: la_dlantp_sub
     public :: la_dlantr
     public :: la_dlantr_sub
#ifdef LA_WITH_XDP
     public :: la_xlangb
     public :: la_xlangb_sub
     public :: la_xlange
     public :: la_xlange_sub
     public :: la_xlangt
     public :: la_xlanhs
     public :: la_xlanhs_sub
     public :: la_xlansb
     public :: la_xlansb_sub
     public :: la_xlansf
     public :: la_xlansf_sub
     public :: la_xlansp
     public :: la_xlansp_sub
     public :: la_xlanst
     public :: la_xlansy
     public :: la_xlansy_sub
     public :: la_xlantb
     public :: la_xlantb_sub
     public :: la_xlantp
     public :: la_xlantp_sub
     public :: la_xlantr
     public :: la_xlantr_sub
#endif
#ifdef LA_WITH_QP
     public :: la_qlangb
     public :: la_qlangb_sub
     public :: la_qlange
     public :: la_qlange_sub
     public :: la_qlangt
     public :: la_qlanhs
     public :: la_qlanhs_sub
     public :: la_qlansb
     public :: la_qlansb_sub
     public :: la_qlansf
     public :: la_qlansf_sub
     public :: la_qlansp
     public :: la_qlansp_sub
     public :: la_qlanst
     public :: la_qlansy
     public :: la_qlansy_sub
     public :: la_qlantb
     public :: la_qlantb_sub
     public :: la_qlantp
     public :: la_qlantp_sub
     public :: la_qlantr
     public :: la_qlantr_sub
#endif
     public :: la_clangb
     public :: la_clangb_sub
     public :: la_clange
     public :: la_clange_sub
     public :: la_clangt
     public :: la_clanhb
     public :: la_clanhb_sub
     public :: la_clanhe
     public :: la_clanhe_sub
     public :: la_clanhf
     public :: la_clanhf_sub
     public :: la_clanhp
     public :: la_clanhp_sub
     public :: la_clanhs
     public :: la_clanhs_sub
     public :: la_clanht
     public :: la_clansb
     public :: la_clansb_sub
     public :: la_clansp
     public :: la_clansp_sub
     public :: la_clansy
     public :: la_clansy_sub
     public :: la_clantb
     public :: la_clantb_sub
     public :: la_clantp
     public :: la_clantp_sub
     public :: la_clantr
     public :: la_clantr_sub
     public :: la_zlangb
     public :: la_zlangb_sub
     public :: la_zlange
     public :: la_zlange_sub
     public :: la_zlangt
     public :: la_zlanhb
     public :: la_zlanhb_sub
     public :: la_zlanhe
     public :: la_zlanhe_sub
     public :: la_zlanhf
     public :: la_zlanhf_sub
     public :: la_zlanhp
     public :: la_zlanhp_sub
     public :: la_zlanhs
     public :: la_zlanhs_sub
     public :: la_zlanht
     public :: la_zlansb
     public :: la_zlansb_sub
     public :: la_zlansp
     public :: la_zlansp_sub
     public :: la_zlansy
     public :: la_zlansy_sub
     public :: la_zlantb
     public :: la_zlantb_sub
     public :: la_zlantp
     public :: la_zlantp_sub
     public :: la_zlantr
     public :: la_zlantr_sub
#ifdef LA_WITH_XDP
     public :: la_ylangb
     public :: la_ylangb_sub
     public :: la_ylange
     public :: la_ylange_sub
     public :: la_ylangt
     public :: la_ylanhb
     public :: la_ylanhb_sub
     public :: la_ylanhe
     public :: la_ylanhe_sub
     public :: la_ylanhf
     public :: la_ylanhf_sub
     public :: la_ylanhp
     public :: la_ylanhp_sub
     public :: la_ylanhs
     public :: la_ylanhs_sub
     public :: la_ylanht
     public :: la_ylansb
     public :: la_ylansb_sub
     public :: la_ylansp
     public :: la_ylansp_sub
     public :: la_ylansy
     public :: la_ylansy_sub
     public :: la_ylantb
     public :: la_ylantb_sub
     public :: la_ylantp
     public :: la_ylantp_sub
     public :: la_ylantr
     public :: la_ylantr_sub
#endif
#ifdef LA_WITH_QP
     public :: la_wlangb
     public :: la_wlangb_sub
     public :: la_wlange
     public :: la_wlange_sub
     public :: la_wlangt
     public :: la_wlanhb
     public :: la_wlanhb_sub
     public :: la_wlanhe
     public :: la_wlanhe_sub
     public :: la_wlanhf
     public :: la_wlanhf_sub
     public :: la_wlanhp
     public :: la_wlanhp_sub
     public :: la_wlanhs
     public :: la_wlanhs_sub
     public :: la_wlanht
     public :: la_wlansb
     public :: la_wlansb_sub
     public :: la_wlansp
     public :: la_wlansp_sub
     public :: la_wlansy
     public :: la_wlansy_sub
     public :: la_wlantb
     public :: la_wlantb_sub
     public :: la_wlantp
     public :: la_wlantp_sub
     public :: la_wlantr
     public :: la_wlantr_sub
#endif

     contains

     !> SLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(sp) function la_slangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           call la_slangb_sub(norm,n,kl,ku,ab,ldab,work,la_slangb)
     end function la_slangb

     !> Pure subroutine form of SLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(sp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_sisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_sisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_slassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slangb_sub
     !> DLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(dp) function la_dlangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           call la_dlangb_sub(norm,n,kl,ku,ab,ldab,work,la_dlangb)
     end function la_dlangb

     !> Pure subroutine form of DLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(dp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_disnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_disnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_dlassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlangb_sub
#ifdef LA_WITH_XDP
     !> XLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(xdp) function la_xlangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           call la_xlangb_sub(norm,n,kl,ku,ab,ldab,work,la_xlangb)
     end function la_xlangb

     !> Pure subroutine form of XLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(xdp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_xisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_xisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_xlassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlangb_sub
#endif
#ifdef LA_WITH_QP
     !> QLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(qp) function la_qlangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           call la_qlangb_sub(norm,n,kl,ku,ab,ldab,work,la_qlangb)
     end function la_qlangb

     !> Pure subroutine form of QLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(qp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_qisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_qisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_qlassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlangb_sub
#endif

     !> SLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(sp) function la_slange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           call la_slange_sub(norm,m,n,a,lda,work,la_slange)
     end function la_slange

     !> Pure subroutine form of SLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slange_sub(norm,m,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_sisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_sisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_slassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slange_sub
     !> DLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(dp) function la_dlange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           call la_dlange_sub(norm,m,n,a,lda,work,la_dlange)
     end function la_dlange

     !> Pure subroutine form of DLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlange_sub(norm,m,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_disnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_disnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_dlassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlange_sub
#ifdef LA_WITH_XDP
     !> XLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(xdp) function la_xlange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           call la_xlange_sub(norm,m,n,a,lda,work,la_xlange)
     end function la_xlange

     !> Pure subroutine form of XLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlange_sub(norm,m,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_xisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_xisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_xlassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlange_sub
#endif
#ifdef LA_WITH_QP
     !> QLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(qp) function la_qlange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           call la_qlange_sub(norm,m,n,a,lda,work,la_qlange)
     end function la_qlange

     !> Pure subroutine form of QLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlange_sub(norm,m,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_qisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_qisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_qlassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlange_sub
#endif

     !> SLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real tridiagonal matrix A.

     pure real(sp) function la_slangt(norm,n,dl,d,du)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_sisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_sisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_sisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_slassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_slassq(n - 1,dl,1,scale,sum)
                 call la_slassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_slangt = anorm
           return
     end function la_slangt
     !> DLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real tridiagonal matrix A.

     pure real(dp) function la_dlangt(norm,n,dl,d,du)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_disnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_disnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_disnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_dlassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_dlassq(n - 1,dl,1,scale,sum)
                 call la_dlassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_dlangt = anorm
           return
     end function la_dlangt
#ifdef LA_WITH_XDP
     !> XLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real tridiagonal matrix A.

     pure real(xdp) function la_xlangt(norm,n,dl,d,du)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_xisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_xisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_xisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_xlassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_xlassq(n - 1,dl,1,scale,sum)
                 call la_xlassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_xlangt = anorm
           return
     end function la_xlangt
#endif
#ifdef LA_WITH_QP
     !> QLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real tridiagonal matrix A.

     pure real(qp) function la_qlangt(norm,n,dl,d,du)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_qisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_qisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_qisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_qlassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_qlassq(n - 1,dl,1,scale,sum)
                 call la_qlassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_qlangt = anorm
           return
     end function la_qlangt
#endif

     !> SLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(sp) function la_slanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           call la_slanhs_sub(norm,n,a,lda,work,la_slanhs)
     end function la_slanhs

     !> Pure subroutine form of SLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slanhs_sub(norm,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_slassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slanhs_sub
     !> DLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(dp) function la_dlanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           call la_dlanhs_sub(norm,n,a,lda,work,la_dlanhs)
     end function la_dlanhs

     !> Pure subroutine form of DLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlanhs_sub(norm,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_dlassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlanhs_sub
#ifdef LA_WITH_XDP
     !> XLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(xdp) function la_xlanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           call la_xlanhs_sub(norm,n,a,lda,work,la_xlanhs)
     end function la_xlanhs

     !> Pure subroutine form of XLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlanhs_sub(norm,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_xlassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlanhs_sub
#endif
#ifdef LA_WITH_QP
     !> QLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(qp) function la_qlanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           call la_qlanhs_sub(norm,n,a,lda,work,la_qlanhs)
     end function la_qlanhs

     !> Pure subroutine form of QLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlanhs_sub(norm,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_qlassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlanhs_sub
#endif

     !> SLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(sp) function la_slansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           call la_slansb_sub(norm,uplo,n,k,ab,ldab,work,la_slansb)
     end function la_slansb

     !> Pure subroutine form of SLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_slassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_slassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_slassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slansb_sub
     !> DLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(dp) function la_dlansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           call la_dlansb_sub(norm,uplo,n,k,ab,ldab,work,la_dlansb)
     end function la_dlansb

     !> Pure subroutine form of DLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_dlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_dlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_dlassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlansb_sub
#ifdef LA_WITH_XDP
     !> XLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(xdp) function la_xlansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           call la_xlansb_sub(norm,uplo,n,k,ab,ldab,work,la_xlansb)
     end function la_xlansb

     !> Pure subroutine form of XLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_xlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_xlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_xlassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlansb_sub
#endif
#ifdef LA_WITH_QP
     !> QLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(qp) function la_qlansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           call la_qlansb_sub(norm,uplo,n,k,ab,ldab,work,la_qlansb)
     end function la_qlansb

     !> Pure subroutine form of QLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_qlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_qlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_qlassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlansb_sub
#endif

     !> SLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(sp) function la_slansf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: a(0:*)
           real(sp),intent(out) :: work(0:*)
           call la_slansf_sub(norm,transr,uplo,n,a,work,la_slansf)
     end function la_slansf

     !> Pure subroutine form of SLANSF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slansf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: a(0:*)
           real(sp),intent(out) :: work(0:*)
           ! Function Value
           real(sp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(sp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(a(0))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='t or 't' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'T')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                 ! a is n by k
                    do j = 0,k - 1
                       do i = 0,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n
                    do j = 0,n - 1
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end do
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    do j = 0,k - 1
                       do i = 0,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n+1
                    do j = 0,n
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              if (ifm == 1) then
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(a(0 + j*lda))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(a(i + j*lda))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(a(0 + j*lda))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                           ! i=j-1-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(a(0))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_slassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_slassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k - 1,a(k),lda + 1,scale,s)
                       ! tri l at a(k,0)
                       call la_slassq(k,a(k - 1),lda + 1,scale,s)
                       ! tri u at a(k-1,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_slassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_slassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(0),lda + 1,scale,s)
                       ! tri l at a(0,0)
                       call la_slassq(k - 1,a(0 + lda),lda + 1,scale,s)
                       ! tri u at a(0,1)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 2
                          call la_slassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_slassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_slassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k - 1,a(0 + k*lda),lda + 1,scale,s)
                       ! tri u at a(0,k)
                       call la_slassq(k,a(0 + (k - 1)*lda),lda + 1,scale,s)
                       ! tri l at a(0,k-1)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_slassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_slassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_slassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                       call la_slassq(k - 1,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_slassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                          ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_slassq(k + j,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(k + 1),lda + 1,scale,s)
                       ! tri l at a(k+1,0)
                       call la_slassq(k,a(k),lda + 1,scale,s)
                       ! tri u at a(k,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_slassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_slassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                       call la_slassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 1
                          call la_slassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                          ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_slassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_slassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                          ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(0 + (k + 1)*lda),lda + 1,scale,s)
                       ! tri u at a(0,k+1)
                       call la_slassq(k,a(0 + k*lda),lda + 1,scale,s)
                       ! tri l at a(0,k)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_slassq(j,a(0 + (j + 1)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_slassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_slassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_slassq(k,a(lda),lda + 1,scale,s)
                       ! tri l at a(0,1)
                       call la_slassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_slansf_sub
     !> DLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(dp) function la_dlansf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: a(0:*)
           real(dp),intent(out) :: work(0:*)
           call la_dlansf_sub(norm,transr,uplo,n,a,work,la_dlansf)
     end function la_dlansf

     !> Pure subroutine form of DLANSF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlansf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: a(0:*)
           real(dp),intent(out) :: work(0:*)
           ! Function Value
           real(dp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(dp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(a(0))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='t or 't' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'T')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                 ! a is n by k
                    do j = 0,k - 1
                       do i = 0,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n
                    do j = 0,n - 1
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end do
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    do j = 0,k - 1
                       do i = 0,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n+1
                    do j = 0,n
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              if (ifm == 1) then
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(a(0 + j*lda))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(a(i + j*lda))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(a(0 + j*lda))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                           ! i=j-1-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(a(0))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_dlassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_dlassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k - 1,a(k),lda + 1,scale,s)
                       ! tri l at a(k,0)
                       call la_dlassq(k,a(k - 1),lda + 1,scale,s)
                       ! tri u at a(k-1,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_dlassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_dlassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(0),lda + 1,scale,s)
                       ! tri l at a(0,0)
                       call la_dlassq(k - 1,a(0 + lda),lda + 1,scale,s)
                       ! tri u at a(0,1)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 2
                          call la_dlassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_dlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_dlassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k - 1,a(0 + k*lda),lda + 1,scale,s)
                       ! tri u at a(0,k)
                       call la_dlassq(k,a(0 + (k - 1)*lda),lda + 1,scale,s)
                       ! tri l at a(0,k-1)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_dlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_dlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_dlassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                       call la_dlassq(k - 1,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_dlassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                          ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_dlassq(k + j,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(k + 1),lda + 1,scale,s)
                       ! tri l at a(k+1,0)
                       call la_dlassq(k,a(k),lda + 1,scale,s)
                       ! tri u at a(k,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_dlassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_dlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                       call la_dlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 1
                          call la_dlassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                          ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_dlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_dlassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                          ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(0 + (k + 1)*lda),lda + 1,scale,s)
                       ! tri u at a(0,k+1)
                       call la_dlassq(k,a(0 + k*lda),lda + 1,scale,s)
                       ! tri l at a(0,k)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_dlassq(j,a(0 + (j + 1)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_dlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_dlassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_dlassq(k,a(lda),lda + 1,scale,s)
                       ! tri l at a(0,1)
                       call la_dlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_dlansf_sub
#ifdef LA_WITH_XDP
     !> XLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(xdp) function la_xlansf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: a(0:*)
           real(xdp),intent(out) :: work(0:*)
           call la_xlansf_sub(norm,transr,uplo,n,a,work,la_xlansf)
     end function la_xlansf

     !> Pure subroutine form of XLANSF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlansf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: a(0:*)
           real(xdp),intent(out) :: work(0:*)
           ! Function Value
           real(xdp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(xdp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(a(0))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='t or 't' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'T')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                 ! a is n by k
                    do j = 0,k - 1
                       do i = 0,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n
                    do j = 0,n - 1
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end do
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    do j = 0,k - 1
                       do i = 0,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n+1
                    do j = 0,n
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              if (ifm == 1) then
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(a(0 + j*lda))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(a(i + j*lda))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(a(0 + j*lda))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                           ! i=j-1-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(a(0))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_xlassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_xlassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k - 1,a(k),lda + 1,scale,s)
                       ! tri l at a(k,0)
                       call la_xlassq(k,a(k - 1),lda + 1,scale,s)
                       ! tri u at a(k-1,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_xlassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_xlassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(0),lda + 1,scale,s)
                       ! tri l at a(0,0)
                       call la_xlassq(k - 1,a(0 + lda),lda + 1,scale,s)
                       ! tri u at a(0,1)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 2
                          call la_xlassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_xlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_xlassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k - 1,a(0 + k*lda),lda + 1,scale,s)
                       ! tri u at a(0,k)
                       call la_xlassq(k,a(0 + (k - 1)*lda),lda + 1,scale,s)
                       ! tri l at a(0,k-1)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_xlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_xlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_xlassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                       call la_xlassq(k - 1,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_xlassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                          ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_xlassq(k + j,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(k + 1),lda + 1,scale,s)
                       ! tri l at a(k+1,0)
                       call la_xlassq(k,a(k),lda + 1,scale,s)
                       ! tri u at a(k,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_xlassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_xlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                       call la_xlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 1
                          call la_xlassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                          ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_xlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_xlassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                          ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(0 + (k + 1)*lda),lda + 1,scale,s)
                       ! tri u at a(0,k+1)
                       call la_xlassq(k,a(0 + k*lda),lda + 1,scale,s)
                       ! tri l at a(0,k)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_xlassq(j,a(0 + (j + 1)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_xlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_xlassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_xlassq(k,a(lda),lda + 1,scale,s)
                       ! tri l at a(0,1)
                       call la_xlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_xlansf_sub
#endif
#ifdef LA_WITH_QP
     !> QLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(qp) function la_qlansf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: a(0:*)
           real(qp),intent(out) :: work(0:*)
           call la_qlansf_sub(norm,transr,uplo,n,a,work,la_qlansf)
     end function la_qlansf

     !> Pure subroutine form of QLANSF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlansf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: a(0:*)
           real(qp),intent(out) :: work(0:*)
           ! Function Value
           real(qp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(qp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(a(0))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='t or 't' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'T')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                 ! a is n by k
                    do j = 0,k - 1
                       do i = 0,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n
                    do j = 0,n - 1
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end do
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    do j = 0,k - 1
                       do i = 0,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end do
                 else
                    ! xpose case; a is k by n+1
                    do j = 0,n
                       do i = 0,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              if (ifm == 1) then
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(a(0 + j*lda))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(a(i + j*lda))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(a(0 + j*lda))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                           ! i=j-1-k
                          aa = abs(a(i + j*lda))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(a(i + j*lda))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(a(0))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(a(i + j*lda))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(a(i + j*lda))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(a(i + j*lda))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_qlassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_qlassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k - 1,a(k),lda + 1,scale,s)
                       ! tri l at a(k,0)
                       call la_qlassq(k,a(k - 1),lda + 1,scale,s)
                       ! tri u at a(k-1,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_qlassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_qlassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(0),lda + 1,scale,s)
                       ! tri l at a(0,0)
                       call la_qlassq(k - 1,a(0 + lda),lda + 1,scale,s)
                       ! tri u at a(0,1)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 2
                          call la_qlassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_qlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_qlassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k - 1,a(0 + k*lda),lda + 1,scale,s)
                       ! tri u at a(0,k)
                       call la_qlassq(k,a(0 + (k - 1)*lda),lda + 1,scale,s)
                       ! tri l at a(0,k-1)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_qlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_qlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_qlassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                       call la_qlassq(k - 1,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_qlassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                          ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_qlassq(k + j,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(k + 1),lda + 1,scale,s)
                       ! tri l at a(k+1,0)
                       call la_qlassq(k,a(k),lda + 1,scale,s)
                       ! tri u at a(k,0)
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_qlassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_qlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(1),lda + 1,scale,s)
                       ! tri l at a(1,0)
                       call la_qlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**t is upper
                       do j = 1,k - 1
                          call la_qlassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                          ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_qlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_qlassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                          ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(0 + (k + 1)*lda),lda + 1,scale,s)
                       ! tri u at a(0,k+1)
                       call la_qlassq(k,a(0 + k*lda),lda + 1,scale,s)
                       ! tri l at a(0,k)
                    else
                       ! a**t is lower
                       do j = 1,k - 1
                          call la_qlassq(j,a(0 + (j + 1)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_qlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_qlassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       call la_qlassq(k,a(lda),lda + 1,scale,s)
                       ! tri l at a(0,1)
                       call la_qlassq(k,a(0),lda + 1,scale,s)
                       ! tri u at a(0,0)
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_qlansf_sub
#endif

     !> SLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(sp) function la_slansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
           call la_slansp_sub(norm,uplo,n,ap,work,la_slansp)
     end function la_slansp

     !> Pure subroutine form of SLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_slassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_slassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (ap(k) /= zero) then
                    absa = abs(ap(k))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slansp_sub
     !> DLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(dp) function la_dlansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
           call la_dlansp_sub(norm,uplo,n,ap,work,la_dlansp)
     end function la_dlansp

     !> Pure subroutine form of DLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_dlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_dlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (ap(k) /= zero) then
                    absa = abs(ap(k))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlansp_sub
#ifdef LA_WITH_XDP
     !> XLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(xdp) function la_xlansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: ap(*)
           real(xdp),intent(out) :: work(*)
           call la_xlansp_sub(norm,uplo,n,ap,work,la_xlansp)
     end function la_xlansp

     !> Pure subroutine form of XLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: ap(*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_xlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_xlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (ap(k) /= zero) then
                    absa = abs(ap(k))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlansp_sub
#endif
#ifdef LA_WITH_QP
     !> QLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(qp) function la_qlansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
           call la_qlansp_sub(norm,uplo,n,ap,work,la_qlansp)
     end function la_qlansp

     !> Pure subroutine form of QLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_qlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_qlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (ap(k) /= zero) then
                    absa = abs(ap(k))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlansp_sub
#endif

     !> SLANST:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric tridiagonal matrix A.

     pure real(sp) function la_slanst(norm,n,d,e)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_slassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_slassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_slanst = anorm
           return
     end function la_slanst
     !> DLANST:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric tridiagonal matrix A.

     pure real(dp) function la_dlanst(norm,n,d,e)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_dlassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_dlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_dlanst = anorm
           return
     end function la_dlanst
#ifdef LA_WITH_XDP
     !> XLANST:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric tridiagonal matrix A.

     pure real(xdp) function la_xlanst(norm,n,d,e)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_xlassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_xlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_xlanst = anorm
           return
     end function la_xlanst
#endif
#ifdef LA_WITH_QP
     !> QLANST:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric tridiagonal matrix A.

     pure real(qp) function la_qlanst(norm,n,d,e)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_qlassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_qlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_qlanst = anorm
           return
     end function la_qlanst
#endif

     !> SLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(sp) function la_slansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           call la_slansy_sub(norm,uplo,n,a,lda,work,la_slansy)
     end function la_slansy

     !> Pure subroutine form of SLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_slassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_slassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_slassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slansy_sub
     !> DLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(dp) function la_dlansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           call la_dlansy_sub(norm,uplo,n,a,lda,work,la_dlansy)
     end function la_dlansy

     !> Pure subroutine form of DLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_dlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_dlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_dlassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlansy_sub
#ifdef LA_WITH_XDP
     !> XLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(xdp) function la_xlansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           call la_xlansy_sub(norm,uplo,n,a,lda,work,la_xlansy)
     end function la_xlansy

     !> Pure subroutine form of XLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_xlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_xlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_xlassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlansy_sub
#endif
#ifdef LA_WITH_QP
     !> QLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(qp) function la_qlansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           call la_qlansy_sub(norm,uplo,n,a,lda,work,la_qlansy)
     end function la_qlansy

     !> Pure subroutine form of QLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_qlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_qlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_qlassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlansy_sub
#endif

     !> SLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(sp) function la_slantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           call la_slantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_slantb)
     end function la_slantb

     !> Pure subroutine form of SLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_slassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_slassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_slassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_slassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slantb_sub
     !> DLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(dp) function la_dlantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           call la_dlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_dlantb)
     end function la_dlantb

     !> Pure subroutine form of DLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_dlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_dlassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_dlassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_dlassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlantb_sub
#ifdef LA_WITH_XDP
     !> XLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(xdp) function la_xlantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           call la_xlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_xlantb)
     end function la_xlantb

     !> Pure subroutine form of XLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_xlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_xlassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_xlassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_xlassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlantb_sub
#endif
#ifdef LA_WITH_QP
     !> QLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(qp) function la_qlantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           call la_qlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_qlantb)
     end function la_qlantb

     !> Pure subroutine form of QLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_qlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_qlassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_qlassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_qlassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlantb_sub
#endif

     !> SLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(sp) function la_slantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
           call la_slantp_sub(norm,uplo,diag,n,ap,work,la_slantp)
     end function la_slantp

     !> Pure subroutine form of SLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_slassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_slassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_slassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_slassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slantp_sub
     !> DLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(dp) function la_dlantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
           call la_dlantp_sub(norm,uplo,diag,n,ap,work,la_dlantp)
     end function la_dlantp

     !> Pure subroutine form of DLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_dlassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_dlassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_dlassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_dlassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlantp_sub
#ifdef LA_WITH_XDP
     !> XLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(xdp) function la_xlantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: ap(*)
           real(xdp),intent(out) :: work(*)
           call la_xlantp_sub(norm,uplo,diag,n,ap,work,la_xlantp)
     end function la_xlantp

     !> Pure subroutine form of XLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: ap(*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_xlassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_xlassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_xlassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_xlassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlantp_sub
#endif
#ifdef LA_WITH_QP
     !> QLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(qp) function la_qlantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
           call la_qlantp_sub(norm,uplo,diag,n,ap,work,la_qlantp)
     end function la_qlantp

     !> Pure subroutine form of QLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_qlassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_qlassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_qlassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_qlassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlantp_sub
#endif

     !> SLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(sp) function la_slantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           call la_slantr_sub(norm,uplo,diag,m,n,a,lda,work,la_slantr)
     end function la_slantr

     !> Pure subroutine form of SLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_slantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_slassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_slassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_slassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_slassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_slantr_sub
     !> DLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(dp) function la_dlantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           call la_dlantr_sub(norm,uplo,diag,m,n,a,lda,work,la_dlantr)
     end function la_dlantr

     !> Pure subroutine form of DLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_dlantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_dlassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_dlassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_dlassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_dlassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_dlantr_sub
#ifdef LA_WITH_XDP
     !> XLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(xdp) function la_xlantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           call la_xlantr_sub(norm,uplo,diag,m,n,a,lda,work,la_xlantr)
     end function la_xlantr

     !> Pure subroutine form of XLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_xlantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_xlassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_xlassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_xlassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_xlassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_xlantr_sub
#endif
#ifdef LA_WITH_QP
     !> QLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(qp) function la_qlantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           call la_qlantr_sub(norm,uplo,diag,m,n,a,lda,work,la_qlantr)
     end function la_qlantr

     !> Pure subroutine form of QLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_qlantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_qlassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_qlassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_qlassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_qlassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_qlantr_sub
#endif

     !> CLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(sp) function la_clangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           call la_clangb_sub(norm,n,kl,ku,ab,ldab,work,la_clangb)
     end function la_clangb

     !> Pure subroutine form of CLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(sp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_sisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_sisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_classq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clangb_sub
     !> ZLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(dp) function la_zlangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           call la_zlangb_sub(norm,n,kl,ku,ab,ldab,work,la_zlangb)
     end function la_zlangb

     !> Pure subroutine form of ZLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(dp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_disnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_disnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_zlassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlangb_sub
#ifdef LA_WITH_XDP
     !> YLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(xdp) function la_ylangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           call la_ylangb_sub(norm,n,kl,ku,ab,ldab,work,la_ylangb)
     end function la_ylangb

     !> Pure subroutine form of YLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(xdp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_xisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_xisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_ylassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylangb_sub
#endif
#ifdef LA_WITH_QP
     !> WLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(qp) function la_wlangb(norm,n,kl,ku,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           call la_wlangb_sub(norm,n,kl,ku,ab,ldab,work,la_wlangb)
     end function la_wlangb

     !> Pure subroutine form of WLANGB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlangb_sub(norm,n,kl,ku,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: kl,ku,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k,l
           real(qp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    temp = abs(ab(i,j))
                    if (value < temp .or. la_qisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                    sum = sum + abs(ab(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 k = ku + 1 - j
                 do i = max(1,j - ku),min(n,j + kl)
                    work(i) = work(i) + abs(ab(k + i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 temp = work(i)
                 if (value < temp .or. la_qisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 l = max(1,j - ku)
                 k = ku + 1 - j + l
                 call la_wlassq(min(n,j + kl) - l + 1,ab(k,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlangb_sub
#endif

     !> CLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(sp) function la_clange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           call la_clange_sub(norm,m,n,a,lda,work,la_clange)
     end function la_clange

     !> Pure subroutine form of CLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clange_sub(norm,m,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_sisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_sisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_classq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clange_sub
     !> ZLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(dp) function la_zlange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           call la_zlange_sub(norm,m,n,a,lda,work,la_zlange)
     end function la_zlange

     !> Pure subroutine form of ZLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlange_sub(norm,m,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_disnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_disnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_zlassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlange_sub
#ifdef LA_WITH_XDP
     !> YLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(xdp) function la_ylange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           call la_ylange_sub(norm,m,n,a,lda,work,la_ylange)
     end function la_ylange

     !> Pure subroutine form of YLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylange_sub(norm,m,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_xisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_xisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_ylassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylange_sub
#endif
#ifdef LA_WITH_QP
     !> WLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(qp) function la_wlange(norm,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           call la_wlange_sub(norm,m,n,a,lda,work,la_wlange)
     end function la_wlange

     !> Pure subroutine form of WLANGE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlange_sub(norm,m,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value,temp
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,m
                    temp = abs(a(i,j))
                    if (value < temp .or. la_qisnan(temp)) value = temp
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,m
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,m
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,m
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,m
                 temp = work(i)
                 if (value < temp .or. la_qisnan(temp)) value = temp
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_wlassq(m,a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlange_sub
#endif

     !> CLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex tridiagonal matrix A.

     pure real(sp) function la_clangt(norm,n,dl,d,du)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_sisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_sisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_sisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_sisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_classq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_classq(n - 1,dl,1,scale,sum)
                 call la_classq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_clangt = anorm
           return
     end function la_clangt
     !> ZLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex tridiagonal matrix A.

     pure real(dp) function la_zlangt(norm,n,dl,d,du)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_disnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_disnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_disnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_disnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_zlassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_zlassq(n - 1,dl,1,scale,sum)
                 call la_zlassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_zlangt = anorm
           return
     end function la_zlangt
#ifdef LA_WITH_XDP
     !> YLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex tridiagonal matrix A.

     pure real(xdp) function la_ylangt(norm,n,dl,d,du)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(xdp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_xisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_xisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_xisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_xisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_ylassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_ylassq(n - 1,dl,1,scale,sum)
                 call la_ylassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_ylangt = anorm
           return
     end function la_ylangt
#endif
#ifdef LA_WITH_QP
     !> WLANGT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex tridiagonal matrix A.

     pure real(qp) function la_wlangt(norm,n,dl,d,du)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(in) :: d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: anorm,scale,sum,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 if (anorm < abs(dl(i)) .or. la_qisnan(abs(dl(i)))) anorm = abs(dl(i))

                 if (anorm < abs(d(i)) .or. la_qisnan(abs(d(i)))) anorm = abs(d(i))

                 if (anorm < abs(du(i)) .or. la_qisnan(abs(du(i)))) anorm = abs(du(i))

              end do
           else if (la_lsame(norm,'O') .or. norm == '1') then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(dl(1))
                 temp = abs(d(n)) + abs(du(n - 1))
                 if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(dl(i)) + abs(du(i - 1))
                    if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(du(1))
                 temp = abs(d(n)) + abs(dl(n - 1))
                 if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 do i = 2,n - 1
                    temp = abs(d(i)) + abs(du(i)) + abs(dl(i - 1))
                    if (anorm < temp .or. la_qisnan(temp)) anorm = temp
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              call la_wlassq(n,d,1,scale,sum)
              if (n > 1) then
                 call la_wlassq(n - 1,dl,1,scale,sum)
                 call la_wlassq(n - 1,du,1,scale,sum)
              end if
              anorm = scale*sqrt(sum)
           end if
           la_wlangt = anorm
           return
     end function la_wlangt
#endif

     !> CLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(sp) function la_clanhb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           call la_clanhb_sub(norm,uplo,n,k,ab,ldab,work,la_clanhb)
     end function la_clanhb

     !> Pure subroutine form of CLANHB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clanhb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    sum = abs(real(ab(k + 1,j),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(ab(1,j),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                    do i = 2,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(ab(k + 1,j),KIND=sp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ab(1,j),KIND=sp))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_classq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_classq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              do j = 1,n
                 if (real(ab(l,j),KIND=sp) /= zero) then
                    absa = abs(real(ab(l,j),KIND=sp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clanhb_sub
     !> ZLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(dp) function la_zlanhb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           call la_zlanhb_sub(norm,uplo,n,k,ab,ldab,work,la_zlanhb)
     end function la_zlanhb

     !> Pure subroutine form of ZLANHB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlanhb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    sum = abs(real(ab(k + 1,j),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(ab(1,j),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                    do i = 2,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(ab(k + 1,j),KIND=dp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ab(1,j),KIND=dp))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_zlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_zlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              do j = 1,n
                 if (real(ab(l,j),KIND=dp) /= zero) then
                    absa = abs(real(ab(l,j),KIND=dp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlanhb_sub
#ifdef LA_WITH_XDP
     !> YLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(xdp) function la_ylanhb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           call la_ylanhb_sub(norm,uplo,n,k,ab,ldab,work,la_ylanhb)
     end function la_ylanhb

     !> Pure subroutine form of YLANHB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylanhb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    sum = abs(real(ab(k + 1,j),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(ab(1,j),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                    do i = 2,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(ab(k + 1,j),KIND=xdp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ab(1,j),KIND=xdp))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_ylassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_ylassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              do j = 1,n
                 if (real(ab(l,j),KIND=xdp) /= zero) then
                    absa = abs(real(ab(l,j),KIND=xdp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylanhb_sub
#endif
#ifdef LA_WITH_QP
     !> WLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(qp) function la_wlanhb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           call la_wlanhb_sub(norm,uplo,n,k,ab,ldab,work,la_wlanhb)
     end function la_wlanhb

     !> Pure subroutine form of WLANHB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlanhb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    sum = abs(real(ab(k + 1,j),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(ab(1,j),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                    do i = 2,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(ab(k + 1,j),KIND=qp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ab(1,j),KIND=qp))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_wlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_wlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              do j = 1,n
                 if (real(ab(l,j),KIND=qp) /= zero) then
                    absa = abs(real(ab(l,j),KIND=qp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlanhb_sub
#endif

     !> CLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(sp) function la_clanhe(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           call la_clanhe_sub(norm,uplo,n,a,lda,work,la_clanhe)
     end function la_clanhe

     !> Pure subroutine form of CLANHE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clanhe_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j - 1
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    sum = abs(real(a(j,j),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(a(j,j),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                    do i = j + 1,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(a(j,j),KIND=sp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(a(j,j),KIND=sp))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_classq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_classq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              do i = 1,n
                 if (real(a(i,i),KIND=sp) /= zero) then
                    absa = abs(real(a(i,i),KIND=sp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clanhe_sub
     !> ZLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(dp) function la_zlanhe(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           call la_zlanhe_sub(norm,uplo,n,a,lda,work,la_zlanhe)
     end function la_zlanhe

     !> Pure subroutine form of ZLANHE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlanhe_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j - 1
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    sum = abs(real(a(j,j),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(a(j,j),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                    do i = j + 1,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(a(j,j),KIND=dp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(a(j,j),KIND=dp))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_zlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_zlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              do i = 1,n
                 if (real(a(i,i),KIND=dp) /= zero) then
                    absa = abs(real(a(i,i),KIND=dp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlanhe_sub
#ifdef LA_WITH_XDP
     !> YLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(xdp) function la_ylanhe(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           call la_ylanhe_sub(norm,uplo,n,a,lda,work,la_ylanhe)
     end function la_ylanhe

     !> Pure subroutine form of YLANHE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylanhe_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j - 1
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    sum = abs(real(a(j,j),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(a(j,j),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                    do i = j + 1,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(a(j,j),KIND=xdp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(a(j,j),KIND=xdp))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_ylassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_ylassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              do i = 1,n
                 if (real(a(i,i),KIND=xdp) /= zero) then
                    absa = abs(real(a(i,i),KIND=xdp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylanhe_sub
#endif
#ifdef LA_WITH_QP
     !> WLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(qp) function la_wlanhe(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           call la_wlanhe_sub(norm,uplo,n,a,lda,work,la_wlanhe)
     end function la_wlanhe

     !> Pure subroutine form of WLANHE: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlanhe_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j - 1
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    sum = abs(real(a(j,j),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    sum = abs(real(a(j,j),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                    do i = j + 1,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(real(a(j,j),KIND=qp))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(a(j,j),KIND=qp))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_wlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_wlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              do i = 1,n
                 if (real(a(i,i),KIND=qp) /= zero) then
                    absa = abs(real(a(i,i),KIND=qp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlanhe_sub
#endif

     !> CLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(sp) function la_clanhf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(0:*)
           complex(sp),intent(in) :: a(0:*)
           call la_clanhf_sub(norm,transr,uplo,n,a,work,la_clanhf)
     end function la_clanhf

     !> Pure subroutine form of CLANHF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clanhf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(0:*)
           complex(sp),intent(in) :: a(0:*)
           ! Function Value
           real(sp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(sp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(real(a(0),KIND=sp))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='c' or 'c' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'C')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is n by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(0,0)
                       temp = abs(real(a(j + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do i = 1,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = k + j - 1
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = k + j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n
                    if (ilu == 1) then
                       ! uplo ='l'
                       do j = 0,k - 2
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j + 1
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       do j = k,n - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       do j = k,n - 1
                          do i = 0,j - k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j - k
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j - k + 1
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j - k + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k)
                       temp = abs(real(a(j + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       temp = abs(real(a(j + 1 + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do i = 2,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j + 1
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = k + j
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j+1; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = k + j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       i = n
                       ! -> u(k-1,k-1)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n+1
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k) at a(0,0)
                       temp = abs(real(a(j + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do j = k + 1,n
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       do j = k + 1,n - 1
                          do i = 0,j - k - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                          i = j - k - 1
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          i = j - k
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=sp))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                          do i = j - k + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_sisnan(temp)) value = temp
                          end do
                       end do
                       j = n
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! u(k,k) at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=sp))
                       if (value < temp .or. la_sisnan(temp)) value = temp
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
             ! find normi(a) ( = norm1(a), since a is hermitian).
              if (ifm == 1) then
                 ! a is 'n'
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(real(a(i + j*lda),KIND=sp))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(real(a(0 + j*lda),KIND=sp))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(real(a(i + j*lda),KIND=sp))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=sp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(real(a(0 + j*lda),KIND=sp))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-1-k
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=sp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(real(a(0),KIND=sp))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=sp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=sp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_sisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_classq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_classq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k - 1
                       ! -> u(k,k) at a(k-1,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=sp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       aa = real(a(l),KIND=sp)
                       ! u(n-1,n-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_classq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 1,k - 2
                          call la_classq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       aa = real(a(0),KIND=sp)
                       ! l(0,0) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(k,k) at a(0,1)
                       do i = 1,k - 1
                          aa = real(a(l),KIND=sp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 2
                          call la_classq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_classq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_classq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda - lda
                       ! -> u(k-1,k-1) at a(0,k-1)
                       aa = real(a(l),KIND=sp)
                       ! u(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k)
                       do j = k,n - 1
                          aa = real(a(l),KIND=sp)
                          ! -> u(j-k,j-k)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_classq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_classq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_classq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(0,0) at a(0,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=sp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! l(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k-1 + (k-1)*lda or l(k-1,k-1) at a(k-1,k-1)
                       aa = real(a(l),KIND=sp)
                       ! l(k-1,k-1) at a(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_classq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                       ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_classq(k + j,a(0 + j*lda),1,scale,s)
                       ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k
                       ! -> u(k,k) at a(k,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=sp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_classq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_classq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=sp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 1
                          call la_classq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                       ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_classq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_classq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                       ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda
                       ! -> u(k,k) at a(0,k)
                       aa = real(a(l),KIND=sp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k+1)
                       do j = k + 1,n - 1
                          aa = real(a(l),KIND=sp)
                          ! -> u(j-k-1,j-k-1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l=k-1+n*lda
                       ! -> u(k-1,k-1) at a(k-1,n)
                       aa = real(a(l),KIND=sp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_classq(j,a(0 + (j + 1)*lda),1,scale,s)
                       ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_classq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_classq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                       ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       aa = real(a(l),KIND=sp)
                       ! l(k,k) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(0,0) at a(0,1)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=sp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=sp)
                          ! l(k+i+1,k+i+1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k - 1 + k*lda or l(k-1,k-1) at a(k-1,k)
                       aa = real(a(l),KIND=sp)
                       ! l(k-1,k-1) at a(k-1,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_clanhf_sub
     !> ZLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(dp) function la_zlanhf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(0:*)
           complex(dp),intent(in) :: a(0:*)
           call la_zlanhf_sub(norm,transr,uplo,n,a,work,la_zlanhf)
     end function la_zlanhf

     !> Pure subroutine form of ZLANHF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlanhf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(0:*)
           complex(dp),intent(in) :: a(0:*)
           ! Function Value
           real(dp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(dp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(real(a(0),KIND=dp))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='c' or 'c' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'C')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is n by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(0,0)
                       temp = abs(real(a(j + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do i = 1,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = k + j - 1
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = i + 1
                          ! =k+j; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = k + j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n
                    if (ilu == 1) then
                       ! uplo ='l'
                       do j = 0,k - 2
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j + 1
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       do j = k,n - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       do j = k,n - 1
                          do i = 0,j - k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j - k
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j - k + 1
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j - k + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k)
                       temp = abs(real(a(j + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       temp = abs(real(a(j + 1 + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do i = 2,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j + 1
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = k + j
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = i + 1
                          ! =k+j+1; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = k + j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       i = n
                       ! -> u(k-1,k-1)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n+1
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k) at a(0,0)
                       temp = abs(real(a(j + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do j = k + 1,n
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       do j = k + 1,n - 1
                          do i = 0,j - k - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                          i = j - k - 1
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          i = j - k
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=dp))
                          if (value < temp .or. la_disnan(temp)) value = temp
                          do i = j - k + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_disnan(temp)) value = temp
                          end do
                       end do
                       j = n
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! u(k,k) at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=dp))
                       if (value < temp .or. la_disnan(temp)) value = temp
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
             ! find normi(a) ( = norm1(a), since a is hermitian).
              if (ifm == 1) then
                 ! a is 'n'
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(real(a(i + j*lda),KIND=dp))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(real(a(0 + j*lda),KIND=dp))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(real(a(i + j*lda),KIND=dp))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=dp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(real(a(0 + j*lda),KIND=dp))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-1-k
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=dp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(real(a(0),KIND=dp))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=dp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=dp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_disnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_zlassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_zlassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k - 1
                       ! -> u(k,k) at a(k-1,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=dp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       aa = real(a(l),KIND=dp)
                       ! u(n-1,n-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_zlassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 1,k - 2
                          call la_zlassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       aa = real(a(0),KIND=dp)
                       ! l(0,0) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(k,k) at a(0,1)
                       do i = 1,k - 1
                          aa = real(a(l),KIND=dp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 2
                          call la_zlassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_zlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_zlassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda - lda
                       ! -> u(k-1,k-1) at a(0,k-1)
                       aa = real(a(l),KIND=dp)
                       ! u(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k)
                       do j = k,n - 1
                          aa = real(a(l),KIND=dp)
                          ! -> u(j-k,j-k)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_zlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_zlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_zlassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(0,0) at a(0,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=dp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! l(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k-1 + (k-1)*lda or l(k-1,k-1) at a(k-1,k-1)
                       aa = real(a(l),KIND=dp)
                       ! l(k-1,k-1) at a(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_zlassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                       ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_zlassq(k + j,a(0 + j*lda),1,scale,s)
                       ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k
                       ! -> u(k,k) at a(k,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=dp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_zlassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_zlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=dp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 1
                          call la_zlassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                       ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_zlassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_zlassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                       ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda
                       ! -> u(k,k) at a(0,k)
                       aa = real(a(l),KIND=dp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k+1)
                       do j = k + 1,n - 1
                          aa = real(a(l),KIND=dp)
                          ! -> u(j-k-1,j-k-1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l=k-1+n*lda
                       ! -> u(k-1,k-1) at a(k-1,n)
                       aa = real(a(l),KIND=dp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_zlassq(j,a(0 + (j + 1)*lda),1,scale,s)
                       ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_zlassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_zlassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                       ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       aa = real(a(l),KIND=dp)
                       ! l(k,k) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(0,0) at a(0,1)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=dp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=dp)
                          ! l(k+i+1,k+i+1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k - 1 + k*lda or l(k-1,k-1) at a(k-1,k)
                       aa = real(a(l),KIND=dp)
                       ! l(k-1,k-1) at a(k-1,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_zlanhf_sub
#ifdef LA_WITH_XDP
     !> YLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(xdp) function la_ylanhf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(0:*)
           complex(xdp),intent(in) :: a(0:*)
           call la_ylanhf_sub(norm,transr,uplo,n,a,work,la_ylanhf)
     end function la_ylanhf

     !> Pure subroutine form of YLANHF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylanhf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(0:*)
           complex(xdp),intent(in) :: a(0:*)
           ! Function Value
           real(xdp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(xdp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(real(a(0),KIND=xdp))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='c' or 'c' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'C')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is n by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(0,0)
                       temp = abs(real(a(j + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do i = 1,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = k + j - 1
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = k + j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n
                    if (ilu == 1) then
                       ! uplo ='l'
                       do j = 0,k - 2
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j + 1
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       do j = k,n - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       do j = k,n - 1
                          do i = 0,j - k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j - k
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j - k + 1
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j - k + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k)
                       temp = abs(real(a(j + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       temp = abs(real(a(j + 1 + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do i = 2,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j + 1
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = k + j
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j+1; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = k + j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       i = n
                       ! -> u(k-1,k-1)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n+1
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k) at a(0,0)
                       temp = abs(real(a(j + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do j = k + 1,n
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       do j = k + 1,n - 1
                          do i = 0,j - k - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                          i = j - k - 1
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          i = j - k
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=xdp))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                          do i = j - k + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_xisnan(temp)) value = temp
                          end do
                       end do
                       j = n
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! u(k,k) at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=xdp))
                       if (value < temp .or. la_xisnan(temp)) value = temp
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
             ! find normi(a) ( = norm1(a), since a is hermitian).
              if (ifm == 1) then
                 ! a is 'n'
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(real(a(i + j*lda),KIND=xdp))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(real(a(0 + j*lda),KIND=xdp))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(real(a(i + j*lda),KIND=xdp))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=xdp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(real(a(0 + j*lda),KIND=xdp))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-1-k
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=xdp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(real(a(0),KIND=xdp))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=xdp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=xdp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_xisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_ylassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_ylassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k - 1
                       ! -> u(k,k) at a(k-1,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=xdp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       aa = real(a(l),KIND=xdp)
                       ! u(n-1,n-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_ylassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 1,k - 2
                          call la_ylassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       aa = real(a(0),KIND=xdp)
                       ! l(0,0) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(k,k) at a(0,1)
                       do i = 1,k - 1
                          aa = real(a(l),KIND=xdp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 2
                          call la_ylassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_ylassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_ylassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda - lda
                       ! -> u(k-1,k-1) at a(0,k-1)
                       aa = real(a(l),KIND=xdp)
                       ! u(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k)
                       do j = k,n - 1
                          aa = real(a(l),KIND=xdp)
                          ! -> u(j-k,j-k)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_ylassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_ylassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_ylassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(0,0) at a(0,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=xdp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! l(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k-1 + (k-1)*lda or l(k-1,k-1) at a(k-1,k-1)
                       aa = real(a(l),KIND=xdp)
                       ! l(k-1,k-1) at a(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_ylassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                       ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_ylassq(k + j,a(0 + j*lda),1,scale,s)
                       ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k
                       ! -> u(k,k) at a(k,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=xdp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_ylassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_ylassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=xdp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 1
                          call la_ylassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                       ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_ylassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_ylassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                       ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda
                       ! -> u(k,k) at a(0,k)
                       aa = real(a(l),KIND=xdp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k+1)
                       do j = k + 1,n - 1
                          aa = real(a(l),KIND=xdp)
                          ! -> u(j-k-1,j-k-1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l=k-1+n*lda
                       ! -> u(k-1,k-1) at a(k-1,n)
                       aa = real(a(l),KIND=xdp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_ylassq(j,a(0 + (j + 1)*lda),1,scale,s)
                       ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_ylassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_ylassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                       ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       aa = real(a(l),KIND=xdp)
                       ! l(k,k) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(0,0) at a(0,1)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=xdp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=xdp)
                          ! l(k+i+1,k+i+1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k - 1 + k*lda or l(k-1,k-1) at a(k-1,k)
                       aa = real(a(l),KIND=xdp)
                       ! l(k-1,k-1) at a(k-1,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_ylanhf_sub
#endif
#ifdef LA_WITH_QP
     !> WLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(qp) function la_wlanhf(norm,transr,uplo,n,a,work)
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(0:*)
           complex(qp),intent(in) :: a(0:*)
           call la_wlanhf_sub(norm,transr,uplo,n,a,work,la_wlanhf)
     end function la_wlanhf

     !> Pure subroutine form of WLANHF: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlanhf_sub(norm,transr,uplo,n,a,work,res)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,transr,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(0:*)
           complex(qp),intent(in) :: a(0:*)
           ! Function Value
           real(qp),intent(out) :: res
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(qp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              res = zero
              return
           else if (n == 1) then
              res = abs(real(a(0),KIND=qp))
              return
           end if
           ! set noe = 1 if n is odd. if n is even set noe=0
           noe = 1
           if (mod(n,2) == 0) noe = 0
           ! set ifm = 0 when form='c' or 'c' and 1 otherwise
           ifm = 1
           if (la_lsame(transr,'C')) ifm = 0
           ! set ilu = 0 when uplo='u or 'u' and 1 otherwise
           ilu = 1
           if (la_lsame(uplo,'U')) ilu = 0
           ! set lda = (n+1)/2 when ifm = 0
           ! set lda = n when ifm = 1 and noe = 1
           ! set lda = n+1 when ifm = 1 and noe = 0
           if (ifm == 1) then
              if (noe == 1) then
                 lda = n
              else
                 ! noe=0
                 lda = n + 1
              end if
           else
              ! ifm=0
              lda = (n + 1)/2
           end if
           if (la_lsame(norm,'M')) then
             ! find max(abs(a(i,j))).
              k = (n + 1)/2
              value = zero
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is n by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(0,0)
                       temp = abs(real(a(j + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do i = 1,n - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = k + j - 1
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = k + j + 1,n - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n
                    if (ilu == 1) then
                       ! uplo ='l'
                       do j = 0,k - 2
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j + 1
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       do j = k,n - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       j = k - 1
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       do j = k,n - 1
                          do i = 0,j - k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j - k
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j - k + 1
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j - k + 2,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is n+1 by k
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k)
                       temp = abs(real(a(j + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       temp = abs(real(a(j + 1 + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do i = 2,n
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j
                          ! l(k+j,k+j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j + 1
                          ! -> l(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 2
                          do i = 0,k + j - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = k + j
                          ! -> u(i,i)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = i + 1
                          ! =k+j+1; i -> u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = k + j + 2,n
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       do i = 0,n - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          ! j=k-1
                       end do
                       ! i=n-1 -> u(n-1,n-1)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       i = n
                       ! -> u(k-1,k-1)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                    end if
                 else
                    ! xpose case; a is k by n+1
                    if (ilu == 1) then
                       ! uplo ='l'
                       j = 0
                       ! -> l(k,k) at a(0,0)
                       temp = abs(real(a(j + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       do j = 1,k - 1
                          do i = 0,j - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j - 1
                          ! l(i,i)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j
                          ! l(j+k,j+k)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! -> l(i,i) is at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do j = k + 1,n
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                    else
                       ! uplo = 'u'
                       do j = 0,k - 1
                          do i = 0,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       j = k
                       ! -> u(j,j) is at a(0,j)
                       temp = abs(real(a(0 + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                       do i = 1,k - 1
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       do j = k + 1,n - 1
                          do i = 0,j - k - 2
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                          i = j - k - 1
                          ! -> u(i,i) at a(i,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          i = j - k
                          ! u(j,j)
                          temp = abs(real(a(i + j*lda),KIND=qp))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                          do i = j - k + 1,k - 1
                             temp = abs(a(i + j*lda))
                             if (value < temp .or. la_qisnan(temp)) value = temp
                          end do
                       end do
                       j = n
                       do i = 0,k - 2
                          temp = abs(a(i + j*lda))
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                       i = k - 1
                       ! u(k,k) at a(i,j)
                       temp = abs(real(a(i + j*lda),KIND=qp))
                       if (value < temp .or. la_qisnan(temp)) value = temp
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
             ! find normi(a) ( = norm1(a), since a is hermitian).
              if (ifm == 1) then
                 ! a is 'n'
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          if (i == k + k) go to 10
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       10 continue
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          if (j > 0) then
                             aa = abs(real(a(i + j*lda),KIND=qp))
                             ! -> a(j+k,j+k)
                             s = s + aa
                             work(i + k) = work(i + k) + s
                             ! i=j
                             i = i + 1
                          end if
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = 0,k - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k + j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(i,j+k)
                             s = s + aa
                             work(i) = work(i) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j+k,j+k)
                          work(j + k) = s + aa
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j,j)
                          work(j) = work(j) + aa
                          s = zero
                          do l = j + 1,k - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu = 1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = k - 1,0,-1
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! -> a(j+k,i+k)
                             s = s + aa
                             work(i + k) = work(i + k) + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j+k,j+k)
                          s = s + aa
                          work(i + k) = work(i + k) + s
                          ! i=j
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! -> a(j,j)
                          work(j) = aa
                          s = zero
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! -> a(l,j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 end if
              else
                 ! ifm=0
                 k = n/2
                 if (noe == 1) then
                    ! n is odd
                    if (ilu == 0) then
                       ! uplo = 'u'
                       n1 = k
                       ! n/2
                       k = k + 1
                       ! k is the row size and lda
                       do i = n1,n - 1
                          work(i) = zero
                       end do
                       do j = 0,n1 - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,n1+i)
                             work(i + n1) = work(i + n1) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=n1=k-1 is special
                       s = abs(real(a(0 + j*lda),KIND=qp))
                       ! a(k-1,k-1)
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k-1,i+n1)
                          work(i + n1) = work(i + n1) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k,n - 1
                          s = zero
                          do i = 0,j - k - 1
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-k
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! a(j-k,j-k)
                          s = s + aa
                          work(j - k) = work(j - k) + s
                          i = i + 1
                          s = abs(real(a(i + j*lda),KIND=qp))
                          ! a(j,j)
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       k = k + 1
                       ! k=(n+1)/2 for n odd and ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 2
                          ! process
                          s = zero
                          do i = 0,j - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! i=j so process of a(j,j)
                          s = s + aa
                          work(j) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k-1 is special :process col a(k-1,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=qp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k,n - 1
                          ! process col j of a = a(j,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 else
                    ! n is even
                    if (ilu == 0) then
                       ! uplo = 'u'
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       do j = 0,k - 1
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j,i+k)
                             work(i + k) = work(i + k) + aa
                             s = s + aa
                          end do
                          work(j) = s
                       end do
                       ! j=k
                       aa = abs(real(a(0 + j*lda),KIND=qp))
                       ! a(k,k)
                       s = aa
                       do i = 1,k - 1
                          aa = abs(a(i + j*lda))
                          ! a(k,k+i)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(j) = work(j) + s
                       do j = k + 1,n - 1
                          s = zero
                          do i = 0,j - 2 - k
                             aa = abs(a(i + j*lda))
                             ! a(i,j-k-1)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          ! i=j-1-k
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! a(j-k-1,j-k-1)
                          s = s + aa
                          work(j - k - 1) = work(j - k - 1) + s
                          i = i + 1
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! a(j,j)
                          s = aa
                          do l = j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(j,l)
                             work(l) = work(l) + aa
                             s = s + aa
                          end do
                          work(j) = work(j) + s
                       end do
                       ! j=n
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(i,k-1)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=qp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = work(i) + s
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    else
                       ! ilu=1
                       do i = k,n - 1
                          work(i) = zero
                       end do
                       ! j=0 is special :process col a(k:n-1,k)
                       s = abs(real(a(0),KIND=qp))
                       ! a(k,k)
                       do i = 1,k - 1
                          aa = abs(a(i))
                          ! a(k+i,k)
                          work(i + k) = work(i + k) + aa
                          s = s + aa
                       end do
                       work(k) = work(k) + s
                       do j = 1,k - 1
                          ! process
                          s = zero
                          do i = 0,j - 2
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          ! i=j-1 so process of a(j-1,j-1)
                          s = s + aa
                          work(j - 1) = s
                          ! is initialised here
                          i = i + 1
                          ! i=j process a(j+k,j+k)
                          aa = abs(real(a(i + j*lda),KIND=qp))
                          s = aa
                          do l = k + j + 1,n - 1
                             i = i + 1
                             aa = abs(a(i + j*lda))
                             ! a(l,k+j)
                             s = s + aa
                             work(l) = work(l) + aa
                          end do
                          work(k + j) = work(k + j) + s
                       end do
                       ! j=k is special :process col a(k,0:k-1)
                       s = zero
                       do i = 0,k - 2
                          aa = abs(a(i + j*lda))
                          ! a(k,i)
                          work(i) = work(i) + aa
                          s = s + aa
                       end do
                       ! i=k-1
                       aa = abs(real(a(i + j*lda),KIND=qp))
                       ! a(k-1,k-1)
                       s = s + aa
                       work(i) = s
                       ! done with col j=k+1
                       do j = k + 1,n
                          ! process col j-1 of a = a(j-1,0:k-1)
                          s = zero
                          do i = 0,k - 1
                             aa = abs(a(i + j*lda))
                             ! a(j-1,i)
                             work(i) = work(i) + aa
                             s = s + aa
                          end do
                          work(j - 1) = work(j - 1) + s
                       end do
                       value = work(0)
                       do i = 1,n - 1
                          temp = work(i)
                          if (value < temp .or. la_qisnan(temp)) value = temp
                       end do
                    end if
                 end if
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
             ! find normf(a).
              k = (n + 1)/2
              scale = zero
              s = one
              if (noe == 1) then
                 ! n is odd
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 3
                          call la_wlassq(k - j - 2,a(k + j + 1 + j*lda),1,scale,s)
                          ! l at a(k,0)
                       end do
                       do j = 0,k - 1
                          call la_wlassq(k + j - 1,a(0 + j*lda),1,scale,s)
                          ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k - 1
                       ! -> u(k,k) at a(k-1,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=qp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       aa = real(a(l),KIND=qp)
                       ! u(n-1,n-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_wlassq(n - j - 1,a(j + 1 + j*lda),1,scale,s)
                          ! trap l at a(0,0)
                       end do
                       do j = 1,k - 2
                          call la_wlassq(j,a(0 + (1 + j)*lda),1,scale,s)
                          ! u at a(0,1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       aa = real(a(0),KIND=qp)
                       ! l(0,0) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(k,k) at a(0,1)
                       do i = 1,k - 1
                          aa = real(a(l),KIND=qp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 2
                          call la_wlassq(j,a(0 + (k + j)*lda),1,scale,s)
                          ! u at a(0,k)
                       end do
                       do j = 0,k - 2
                          call la_wlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_wlassq(k - j - 1,a(j + 1 + (j + k - 1)*lda),1,scale,s)
                          ! l at a(0,k-1)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda - lda
                       ! -> u(k-1,k-1) at a(0,k-1)
                       aa = real(a(l),KIND=qp)
                       ! u(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k)
                       do j = k,n - 1
                          aa = real(a(l),KIND=qp)
                          ! -> u(j-k,j-k)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_wlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       do j = k,n - 1
                          call la_wlassq(k,a(0 + j*lda),1,scale,s)
                          ! k by k-1 rect. at a(0,k)
                       end do
                       do j = 0,k - 3
                          call la_wlassq(k - j - 2,a(j + 2 + j*lda),1,scale,s)
                          ! l at a(1,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(0,0) at a(0,0)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=qp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! l(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k-1 + (k-1)*lda or l(k-1,k-1) at a(k-1,k-1)
                       aa = real(a(l),KIND=qp)
                       ! l(k-1,k-1) at a(k-1,k-1)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              else
                 ! n is even
                 if (ifm == 1) then
                    ! a is normal
                    if (ilu == 0) then
                       ! a is upper
                       do j = 0,k - 2
                          call la_wlassq(k - j - 1,a(k + j + 2 + j*lda),1,scale,s)
                       ! l at a(k+1,0)
                       end do
                       do j = 0,k - 1
                          call la_wlassq(k + j,a(0 + j*lda),1,scale,s)
                       ! trap u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = k
                       ! -> u(k,k) at a(k,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=qp)
                          ! u(k+i,k+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! u(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    else
                       ! ilu=1
                       do j = 0,k - 1
                          call la_wlassq(n - j - 1,a(j + 2 + j*lda),1,scale,s)
                          ! trap l at a(1,0)
                       end do
                       do j = 1,k - 1
                          call la_wlassq(j,a(0 + j*lda),1,scale,s)
                          ! u at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       do i = 0,k - 1
                          aa = real(a(l),KIND=qp)
                          ! l(k-1+i,k-1+i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                    end if
                 else
                    ! a is xpose
                    if (ilu == 0) then
                       ! a**h is upper
                       do j = 1,k - 1
                          call la_wlassq(j,a(0 + (k + 1 + j)*lda),1,scale,s)
                       ! u at a(0,k+1)
                       end do
                       do j = 0,k - 1
                          call la_wlassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,0)
                       end do
                       do j = 0,k - 2
                          call la_wlassq(k - j - 1,a(j + 1 + (j + k)*lda),1,scale,s)
                       ! l at a(0,k)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0 + k*lda
                       ! -> u(k,k) at a(0,k)
                       aa = real(a(l),KIND=qp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = l + lda
                       ! -> u(0,0) at a(0,k+1)
                       do j = k + 1,n - 1
                          aa = real(a(l),KIND=qp)
                          ! -> u(j-k-1,j-k-1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! -> u(j,j)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l=k-1+n*lda
                       ! -> u(k-1,k-1) at a(k-1,n)
                       aa = real(a(l),KIND=qp)
                       ! u(k,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    else
                       ! a**h is lower
                       do j = 1,k - 1
                          call la_wlassq(j,a(0 + (j + 1)*lda),1,scale,s)
                       ! u at a(0,1)
                       end do
                       do j = k + 1,n
                          call la_wlassq(k,a(0 + j*lda),1,scale,s)
                       ! k by k rect. at a(0,k+1)
                       end do
                       do j = 0,k - 2
                          call la_wlassq(k - j - 1,a(j + 1 + j*lda),1,scale,s)
                       ! l at a(0,0)
                       end do
                       s = s + s
                       ! double s for the off diagonal elements
                       l = 0
                       ! -> l(k,k) at a(0,0)
                       aa = real(a(l),KIND=qp)
                       ! l(k,k) at a(0,0)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                       l = lda
                       ! -> l(0,0) at a(0,1)
                       do i = 0,k - 2
                          aa = real(a(l),KIND=qp)
                          ! l(i,i)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          aa = real(a(l + 1),KIND=qp)
                          ! l(k+i+1,k+i+1)
                          if (aa /= zero) then
                             if (scale < aa) then
                                s = one + s*(scale/aa)**2
                                scale = aa
                             else
                                s = s + (aa/scale)**2
                             end if
                          end if
                          l = l + lda + 1
                       end do
                       ! l-> k - 1 + k*lda or l(k-1,k-1) at a(k-1,k)
                       aa = real(a(l),KIND=qp)
                       ! l(k-1,k-1) at a(k-1,k)
                       if (aa /= zero) then
                          if (scale < aa) then
                             s = one + s*(scale/aa)**2
                             scale = aa
                          else
                             s = s + (aa/scale)**2
                          end if
                       end if
                    end if
                 end if
              end if
              value = scale*sqrt(s)
           end if
           res = value
           return
     end subroutine la_wlanhf_sub
#endif

     !> CLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(sp) function la_clanhp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           call la_clanhp_sub(norm,uplo,n,ap,work,la_clanhp)
     end function la_clanhp

     !> Pure subroutine form of CLANHP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clanhp_sub(norm,uplo,n,ap,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 0
                 do j = 1,n
                    do i = k + 1,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + j
                    sum = abs(real(ap(k),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 k = 1
                 do j = 1,n
                    sum = abs(real(ap(k),KIND=sp))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                    do i = k + 1,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(real(ap(k),KIND=sp))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ap(k),KIND=sp))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_classq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_classq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=sp) /= zero) then
                    absa = abs(real(ap(k),KIND=sp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clanhp_sub
     !> ZLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(dp) function la_zlanhp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           call la_zlanhp_sub(norm,uplo,n,ap,work,la_zlanhp)
     end function la_zlanhp

     !> Pure subroutine form of ZLANHP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlanhp_sub(norm,uplo,n,ap,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 0
                 do j = 1,n
                    do i = k + 1,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + j
                    sum = abs(real(ap(k),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 k = 1
                 do j = 1,n
                    sum = abs(real(ap(k),KIND=dp))
                    if (value < sum .or. la_disnan(sum)) value = sum
                    do i = k + 1,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(real(ap(k),KIND=dp))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ap(k),KIND=dp))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_zlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_zlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=dp) /= zero) then
                    absa = abs(real(ap(k),KIND=dp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlanhp_sub
#ifdef LA_WITH_XDP
     !> YLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(xdp) function la_ylanhp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           call la_ylanhp_sub(norm,uplo,n,ap,work,la_ylanhp)
     end function la_ylanhp

     !> Pure subroutine form of YLANHP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylanhp_sub(norm,uplo,n,ap,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 0
                 do j = 1,n
                    do i = k + 1,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + j
                    sum = abs(real(ap(k),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 k = 1
                 do j = 1,n
                    sum = abs(real(ap(k),KIND=xdp))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                    do i = k + 1,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(real(ap(k),KIND=xdp))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ap(k),KIND=xdp))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_ylassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_ylassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=xdp) /= zero) then
                    absa = abs(real(ap(k),KIND=xdp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylanhp_sub
#endif
#ifdef LA_WITH_QP
     !> WLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(qp) function la_wlanhp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           call la_wlanhp_sub(norm,uplo,n,ap,work,la_wlanhp)
     end function la_wlanhp

     !> Pure subroutine form of WLANHP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlanhp_sub(norm,uplo,n,ap,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 0
                 do j = 1,n
                    do i = k + 1,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + j
                    sum = abs(real(ap(k),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 k = 1
                 do j = 1,n
                    sum = abs(real(ap(k),KIND=qp))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                    do i = k + 1,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is hermitian).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(real(ap(k),KIND=qp))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(real(ap(k),KIND=qp))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_wlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_wlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=qp) /= zero) then
                    absa = abs(real(ap(k),KIND=qp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlanhp_sub
#endif

     !> CLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(sp) function la_clanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           call la_clanhs_sub(norm,n,a,lda,work,la_clanhs)
     end function la_clanhs

     !> Pure subroutine form of CLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clanhs_sub(norm,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_classq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clanhs_sub
     !> ZLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(dp) function la_zlanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           call la_zlanhs_sub(norm,n,a,lda,work,la_zlanhs)
     end function la_zlanhs

     !> Pure subroutine form of ZLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlanhs_sub(norm,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_zlassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlanhs_sub
#ifdef LA_WITH_XDP
     !> YLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(xdp) function la_ylanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           call la_ylanhs_sub(norm,n,a,lda,work,la_ylanhs)
     end function la_ylanhs

     !> Pure subroutine form of YLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylanhs_sub(norm,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_ylassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylanhs_sub
#endif
#ifdef LA_WITH_QP
     !> WLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(qp) function la_wlanhs(norm,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           call la_wlanhs_sub(norm,n,a,lda,work,la_wlanhs)
     end function la_wlanhs

     !> Pure subroutine form of WLANHS: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlanhs_sub(norm,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    sum = abs(a(i,j))
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end do
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              do j = 1,n
                 sum = zero
                 do i = 1,min(n,j + 1)
                    sum = sum + abs(a(i,j))
                 end do
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              do i = 1,n
                 work(i) = zero
              end do
              do j = 1,n
                 do i = 1,min(n,j + 1)
                    work(i) = work(i) + abs(a(i,j))
                 end do
              end do
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              do j = 1,n
                 call la_wlassq(min(n,j + 1),a(1,j),1,scale,sum)
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlanhs_sub
#endif

     !> CLANHT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian tridiagonal matrix A.

     pure real(sp) function la_clanht(norm,n,d,e)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: d(*)
           complex(sp),intent(in) :: e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_sisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_classq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_slassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_clanht = anorm
           return
     end function la_clanht
     !> ZLANHT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian tridiagonal matrix A.

     pure real(dp) function la_zlanht(norm,n,d,e)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: d(*)
           complex(dp),intent(in) :: e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_disnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_zlassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_dlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_zlanht = anorm
           return
     end function la_zlanht
#ifdef LA_WITH_XDP
     !> YLANHT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian tridiagonal matrix A.

     pure real(xdp) function la_ylanht(norm,n,d,e)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(in) :: d(*)
           complex(xdp),intent(in) :: e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(xdp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_xisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_ylassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_xlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_ylanht = anorm
           return
     end function la_ylanht
#endif
#ifdef LA_WITH_QP
     !> WLANHT:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian tridiagonal matrix A.

     pure real(qp) function la_wlanht(norm,n,d,e)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: d(*)
           complex(qp),intent(in) :: e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: anorm,scale,sum
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n <= 0) then
              anorm = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              anorm = abs(d(n))
              do i = 1,n - 1
                 sum = abs(d(i))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 sum = abs(e(i))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
              end do
           else if (la_lsame(norm,'O') .or. norm == '1' .or. la_lsame(norm,'I')) &
                     then
              ! find norm1(a).
              if (n == 1) then
                 anorm = abs(d(1))
              else
                 anorm = abs(d(1)) + abs(e(1))
                 sum = abs(e(n - 1)) + abs(d(n))
                 if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 do i = 2,n - 1
                    sum = abs(d(i)) + abs(e(i)) + abs(e(i - 1))
                    if (anorm < sum .or. la_qisnan(sum)) anorm = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (n > 1) then
                 call la_wlassq(n - 1,e,1,scale,sum)
                 sum = 2*sum
              end if
              call la_qlassq(n,d,1,scale,sum)
              anorm = scale*sqrt(sum)
           end if
           la_wlanht = anorm
           return
     end function la_wlanht
#endif

     !> CLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(sp) function la_clansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           call la_clansb_sub(norm,uplo,n,k,ab,ldab,work,la_clansb)
     end function la_clansb

     !> Pure subroutine form of CLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_classq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_classq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_classq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clansb_sub
     !> ZLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(dp) function la_zlansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           call la_zlansb_sub(norm,uplo,n,k,ab,ldab,work,la_zlansb)
     end function la_zlansb

     !> Pure subroutine form of ZLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_zlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_zlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_zlassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlansb_sub
#ifdef LA_WITH_XDP
     !> YLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(xdp) function la_ylansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           call la_ylansb_sub(norm,uplo,n,k,ab,ldab,work,la_ylansb)
     end function la_ylansb

     !> Pure subroutine form of YLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_ylassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_ylassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_ylassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylansb_sub
#endif
#ifdef LA_WITH_QP
     !> WLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(qp) function la_wlansb(norm,uplo,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           call la_wlansb_sub(norm,uplo,n,k,ab,ldab,work,la_wlansb)
     end function la_wlansb

     !> Pure subroutine form of WLANSB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlansb_sub(norm,uplo,n,k,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = max(k + 2 - j,1),k + 1
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = 1,min(n + 1 - j,k + 1)
                       sum = abs(ab(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    l = k + 1 - j
                    do i = max(1,j - k),j - 1
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(ab(k + 1,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ab(1,j))
                    l = 1 - j
                    do i = j + 1,min(n,j + k)
                       absa = abs(ab(l + i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (k > 0) then
                 if (la_lsame(uplo,'U')) then
                    do j = 2,n
                       call la_wlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                    l = k + 1
                 else
                    do j = 1,n - 1
                       call la_wlassq(min(n - j,k),ab(2,j),1,scale,sum)
                    end do
                    l = 1
                 end if
                 sum = 2*sum
              else
                 l = 1
              end if
              call la_wlassq(n,ab(l,1),ldab,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlansb_sub
#endif

     !> CLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(sp) function la_clansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           call la_clansp_sub(norm,uplo,n,ap,work,la_clansp)
     end function la_clansp

     !> Pure subroutine form of CLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,aimag,real,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_classq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_classq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=sp) /= zero) then
                    absa = abs(real(ap(k),KIND=sp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (aimag(ap(k)) /= zero) then
                    absa = abs(aimag(ap(k)))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clansp_sub
     !> ZLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(dp) function la_zlansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           call la_zlansp_sub(norm,uplo,n,ap,work,la_zlansp)
     end function la_zlansp

     !> Pure subroutine form of ZLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_zlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_zlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=dp) /= zero) then
                    absa = abs(real(ap(k),KIND=dp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (aimag(ap(k)) /= zero) then
                    absa = abs(aimag(ap(k)))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlansp_sub
#ifdef LA_WITH_XDP
     !> YLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(xdp) function la_ylansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           call la_ylansp_sub(norm,uplo,n,ap,work,la_ylansp)
     end function la_ylansp

     !> Pure subroutine form of YLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_ylassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_ylassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=xdp) /= zero) then
                    absa = abs(real(ap(k),KIND=xdp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (aimag(ap(k)) /= zero) then
                    absa = abs(aimag(ap(k)))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylansp_sub
#endif
#ifdef LA_WITH_QP
     !> WLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(qp) function la_wlansp(norm,uplo,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           call la_wlansp_sub(norm,uplo,n,ap,work,la_wlansp)
     end function la_wlansp

     !> Pure subroutine form of WLANSP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlansp_sub(norm,uplo,n,ap,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,k
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 k = 1
                 do j = 1,n
                    do i = k,k + j - 1
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + j
                 end do
              else
                 k = 1
                 do j = 1,n
                    do i = k,k + n - j
                       sum = abs(ap(i))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                    k = k + n - j + 1
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              k = 1
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    work(j) = sum + abs(ap(k))
                    k = k + 1
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(ap(k))
                    k = k + 1
                    do i = j + 1,n
                       absa = abs(ap(k))
                       sum = sum + absa
                       work(i) = work(i) + absa
                       k = k + 1
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              k = 2
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_wlassq(j - 1,ap(k),1,scale,sum)
                    k = k + j
                 end do
              else
                 do j = 1,n - 1
                    call la_wlassq(n - j,ap(k),1,scale,sum)
                    k = k + n - j + 1
                 end do
              end if
              sum = 2*sum
              k = 1
              do i = 1,n
                 if (real(ap(k),KIND=qp) /= zero) then
                    absa = abs(real(ap(k),KIND=qp))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (aimag(ap(k)) /= zero) then
                    absa = abs(aimag(ap(k)))
                    if (scale < absa) then
                       sum = one + sum*(scale/absa)**2
                       scale = absa
                    else
                       sum = sum + (absa/scale)**2
                    end if
                 end if
                 if (la_lsame(uplo,'U')) then
                    k = k + i + 1
                 else
                    k = k + n - i + 1
                 end if
              end do
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlansp_sub
#endif

     !> CLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(sp) function la_clansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           call la_clansy_sub(norm,uplo,n,a,lda,work,la_clansy)
     end function la_clansy

     !> Pure subroutine form of CLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_sisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_classq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_classq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_classq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clansy_sub
     !> ZLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(dp) function la_zlansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           call la_zlansy_sub(norm,uplo,n,a,lda,work,la_zlansy)
     end function la_zlansy

     !> Pure subroutine form of ZLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_disnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_zlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_zlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_zlassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlansy_sub
#ifdef LA_WITH_XDP
     !> YLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(xdp) function la_ylansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           call la_ylansy_sub(norm,uplo,n,a,lda,work,la_ylansy)
     end function la_ylansy

     !> Pure subroutine form of YLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_xisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_ylassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_ylassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_ylassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylansy_sub
#endif
#ifdef LA_WITH_QP
     !> WLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(qp) function la_wlansy(norm,uplo,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           call la_wlansy_sub(norm,uplo,n,a,lda,work,la_wlansy)
     end function la_wlansy

     !> Pure subroutine form of WLANSY: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlansy_sub(norm,uplo,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: norm,uplo
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: absa,scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    do i = 1,j
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              else
                 do j = 1,n
                    do i = j,n
                       sum = abs(a(i,j))
                       if (value < sum .or. la_qisnan(sum)) value = sum
                    end do
                 end do
              end if
           else if ((la_lsame(norm,'I')) .or. (la_lsame(norm,'O')) .or. ( &
                     norm == '1')) then
              ! find normi(a) ( = norm1(a), since a is symmetric).
              value = zero
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    sum = zero
                    do i = 1,j - 1
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    work(j) = sum + abs(a(j,j))
                 end do
                 do i = 1,n
                    sum = work(i)
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do i = 1,n
                    work(i) = zero
                 end do
                 do j = 1,n
                    sum = work(j) + abs(a(j,j))
                    do i = j + 1,n
                       absa = abs(a(i,j))
                       sum = sum + absa
                       work(i) = work(i) + absa
                    end do
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              scale = zero
              sum = one
              if (la_lsame(uplo,'U')) then
                 do j = 2,n
                    call la_wlassq(j - 1,a(1,j),1,scale,sum)
                 end do
              else
                 do j = 1,n - 1
                    call la_wlassq(n - j,a(j + 1,j),1,scale,sum)
                 end do
              end if
              sum = 2*sum
              call la_wlassq(n,a,lda + 1,scale,sum)
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlansy_sub
#endif

     !> CLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(sp) function la_clantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           call la_clantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_clantb)
     end function la_clantb

     !> Pure subroutine form of CLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_classq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_classq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_classq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_classq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clantb_sub
     !> ZLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(dp) function la_zlantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           call la_zlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_zlantb)
     end function la_zlantb

     !> Pure subroutine form of ZLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_zlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_zlassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_zlassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_zlassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlantb_sub
#ifdef LA_WITH_XDP
     !> YLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(xdp) function la_ylantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           call la_ylantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_ylantb)
     end function la_ylantb

     !> Pure subroutine form of YLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_ylassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_ylassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_ylassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_ylassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylantb_sub
#endif
#ifdef LA_WITH_QP
     !> WLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(qp) function la_wlantb(norm,uplo,diag,n,k,ab,ldab,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           call la_wlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,la_wlantb)
     end function la_wlantb

     !> Pure subroutine form of WLANTB: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlantb_sub(norm,uplo,diag,n,k,ab,ldab,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: k,ldab,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ab(ldab,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,l
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = max(k + 2 - j,1),k + 1
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = abs(ab(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = max(k + 2 - j,1),k
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = max(k + 2 - j,1),k + 1
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = 2,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(n + 1 - j,k + 1)
                          sum = sum + abs(ab(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              value = zero
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j - 1
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = k + 1 - j
                       do i = max(1,j - k),j
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j + 1,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       l = 1 - j
                       do i = j,min(n,j + k)
                          work(i) = work(i) + abs(ab(l + i,j))
                       end do
                    end do
                 end if
              end if
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 2,n
                          call la_wlassq(min(j - 1,k),ab(max(k + 2 - j,1),j),1,scale, &
                                    sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_wlassq(min(j,k + 1),ab(max(k + 2 - j,1),j),1,scale,sum)

                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    if (k > 0) then
                       do j = 1,n - 1
                          call la_wlassq(min(n - j,k),ab(2,j),1,scale,sum)
                       end do
                    end if
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_wlassq(min(n - j + 1,k + 1),ab(1,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlantb_sub
#endif

     !> CLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(sp) function la_clantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           call la_clantp_sub(norm,uplo,diag,n,ap,work,la_clantp)
     end function la_clantp

     !> Pure subroutine form of CLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: ap(*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_classq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_classq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_classq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_classq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clantp_sub
     !> ZLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(dp) function la_zlantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           call la_zlantp_sub(norm,uplo,diag,n,ap,work,la_zlantp)
     end function la_zlantp

     !> Pure subroutine form of ZLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: ap(*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_zlassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_zlassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_zlassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_zlassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlantp_sub
#ifdef LA_WITH_XDP
     !> YLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(xdp) function la_ylantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           call la_ylantp_sub(norm,uplo,diag,n,ap,work,la_ylantp)
     end function la_ylantp

     !> Pure subroutine form of YLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: ap(*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_ylassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_ylassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_ylassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_ylassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylantp_sub
#endif
#ifdef LA_WITH_QP
     !> WLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(qp) function la_wlantp(norm,uplo,diag,n,ap,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           call la_wlantp_sub(norm,uplo,diag,n,ap,work,la_wlantp)
     end function la_wlantp

     !> Pure subroutine form of WLANTP: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlantp_sub(norm,uplo,diag,n,ap,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: ap(*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j,k
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              k = 1
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 2
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k + 1,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = k,k + j - 1
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + j
                    end do
                 else
                    do j = 1,n
                       do i = k,k + n - j
                          sum = abs(ap(i))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                       k = k + n - j + 1
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              k = 1
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k,k + j - 2
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + j - 1
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + j
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = k + 1,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    else
                       sum = zero
                       do i = k,k + n - j
                          sum = sum + abs(ap(i))
                       end do
                    end if
                    k = k + n - j + 1
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              k = 1
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,j - 1
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                       k = k + 1
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,j
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,n
                       work(i) = one
                    end do
                    do j = 1,n
                       k = k + 1
                       do i = j + 1,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 else
                    do i = 1,n
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,n
                          work(i) = work(i) + abs(ap(k))
                          k = k + 1
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,n
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 2,n
                       call la_wlassq(j - 1,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_wlassq(j,ap(k),1,scale,sum)
                       k = k + j
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = n
                    k = 2
                    do j = 1,n - 1
                       call la_wlassq(n - j,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 else
                    scale = zero
                    sum = one
                    k = 1
                    do j = 1,n
                       call la_wlassq(n - j + 1,ap(k),1,scale,sum)
                       k = k + n - j + 1
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlantp_sub
#endif

     !> CLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(sp) function la_clantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           call la_clantr_sub(norm,uplo,diag,m,n,a,lda,work,la_clantr)
     end function la_clantr

     !> Pure subroutine form of CLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_clantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(out) :: work(*)
           complex(sp),intent(in) :: a(lda,*)
           ! Function Value
           real(sp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(sp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_sisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_sisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_sisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_classq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_classq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_classq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_classq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_clantr_sub
     !> ZLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(dp) function la_zlantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           call la_zlantr_sub(norm,uplo,diag,m,n,a,lda,work,la_zlantr)
     end function la_zlantr

     !> Pure subroutine form of ZLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_zlantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(out) :: work(*)
           complex(dp),intent(in) :: a(lda,*)
           ! Function Value
           real(dp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(dp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_disnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_disnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_disnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_zlassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_zlassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_zlassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_zlassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_zlantr_sub
#ifdef LA_WITH_XDP
     !> YLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(xdp) function la_ylantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           call la_ylantr_sub(norm,uplo,diag,m,n,a,lda,work,la_ylantr)
     end function la_ylantr

     !> Pure subroutine form of YLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_ylantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(in) :: a(lda,*)
           ! Function Value
           real(xdp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(xdp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_xisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_xisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_xisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_ylassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_ylassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_ylassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_ylassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_ylantr_sub
#endif
#ifdef LA_WITH_QP
     !> WLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(qp) function la_wlantr(norm,uplo,diag,m,n,a,lda,work)
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           call la_wlantr_sub(norm,uplo,diag,m,n,a,lda,work,la_wlantr)
     end function la_wlantr

     !> Pure subroutine form of WLANTR: the same arguments, followed by RES,
     !> which returns the value of the function.
     pure subroutine la_wlantr_sub(norm,uplo,diag,m,n,a,lda,work,res)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(out) :: work(*)
           complex(qp),intent(in) :: a(lda,*)
           ! Function Value
           real(qp),intent(out) :: res
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: udiag
           integer(ilp) :: i,j
           real(qp) :: scale,sum,value
           ! Intrinsic Functions
           intrinsic :: abs,min,sqrt
           ! Executable Statements
           if (min(m,n) == 0) then
              value = zero
           else if (la_lsame(norm,'M')) then
              ! find max(abs(a(i,j))).
              if (la_lsame(diag,'U')) then
                 value = one
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j + 1,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              else
                 value = zero
                 if (la_lsame(uplo,'U')) then
                    do j = 1,n
                       do i = 1,min(m,j)
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 else
                    do j = 1,n
                       do i = j,m
                          sum = abs(a(i,j))
                          if (value < sum .or. la_qisnan(sum)) value = sum
                       end do
                    end do
                 end if
              end if
           else if ((la_lsame(norm,'O')) .or. (norm == '1')) then
              ! find norm1(a).
              value = zero
              udiag = la_lsame(diag,'U')
              if (la_lsame(uplo,'U')) then
                 do j = 1,n
                    if ((udiag) .and. (j <= m)) then
                       sum = one
                       do i = 1,j - 1
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = 1,min(m,j)
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              else
                 do j = 1,n
                    if (udiag) then
                       sum = one
                       do i = j + 1,m
                          sum = sum + abs(a(i,j))
                       end do
                    else
                       sum = zero
                       do i = j,m
                          sum = sum + abs(a(i,j))
                       end do
                    end if
                    if (value < sum .or. la_qisnan(sum)) value = sum
                 end do
              end if
           else if (la_lsame(norm,'I')) then
              ! find normi(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    do i = 1,m
                       work(i) = one
                    end do
                    do j = 1,n
                       do i = 1,min(m,j - 1)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = 1,min(m,j)
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    do i = 1,min(m,n)
                       work(i) = one
                    end do
                    do i = n + 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j + 1,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 else
                    do i = 1,m
                       work(i) = zero
                    end do
                    do j = 1,n
                       do i = j,m
                          work(i) = work(i) + abs(a(i,j))
                       end do
                    end do
                 end if
              end if
              value = zero
              do i = 1,m
                 sum = work(i)
                 if (value < sum .or. la_qisnan(sum)) value = sum
              end do
           else if ((la_lsame(norm,'F')) .or. (la_lsame(norm,'E'))) &
                     then
              ! find normf(a).
              if (la_lsame(uplo,'U')) then
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 2,n
                       call la_wlassq(min(m,j - 1),a(1,j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_wlassq(min(m,j),a(1,j),1,scale,sum)
                    end do
                 end if
              else
                 if (la_lsame(diag,'U')) then
                    scale = one
                    sum = min(m,n)
                    do j = 1,n
                       call la_wlassq(m - j,a(min(m,j + 1),j),1,scale,sum)
                    end do
                 else
                    scale = zero
                    sum = one
                    do j = 1,n
                       call la_wlassq(m - j + 1,a(j,j),1,scale,sum)
                    end do
                 end if
              end if
              value = scale*sqrt(sum)
           end if
           res = value
           return
     end subroutine la_wlantr_sub
#endif

end module la_lapack_blas_like_mnorm
