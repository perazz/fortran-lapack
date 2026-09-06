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
     public :: la_slange
     public :: la_slangt
     public :: la_slanhs
     public :: la_slansb
     public :: la_slansf
     public :: la_slansp
     public :: la_slanst
     public :: la_slansy
     public :: la_slantb
     public :: la_slantp
     public :: la_slantr
     public :: la_dlangb
     public :: la_dlange
     public :: la_dlangt
     public :: la_dlanhs
     public :: la_dlansb
     public :: la_dlansf
     public :: la_dlansp
     public :: la_dlanst
     public :: la_dlansy
     public :: la_dlantb
     public :: la_dlantp
     public :: la_dlantr
     public :: la_qlangb
     public :: la_qlange
     public :: la_qlangt
     public :: la_qlanhs
     public :: la_qlansb
     public :: la_qlansf
     public :: la_qlansp
     public :: la_qlanst
     public :: la_qlansy
     public :: la_qlantb
     public :: la_qlantp
     public :: la_qlantr
     public :: la_clangb
     public :: la_clange
     public :: la_clangt
     public :: la_clanhb
     public :: la_clanhe
     public :: la_clanhf
     public :: la_clanhp
     public :: la_clanhs
     public :: la_clanht
     public :: la_clansb
     public :: la_clansp
     public :: la_clansy
     public :: la_clantb
     public :: la_clantp
     public :: la_clantr
     public :: la_zlangb
     public :: la_zlange
     public :: la_zlangt
     public :: la_zlanhb
     public :: la_zlanhe
     public :: la_zlanhf
     public :: la_zlanhp
     public :: la_zlanhs
     public :: la_zlanht
     public :: la_zlansb
     public :: la_zlansp
     public :: la_zlansy
     public :: la_zlantb
     public :: la_zlantp
     public :: la_zlantr
     public :: la_wlangb
     public :: la_wlange
     public :: la_wlangt
     public :: la_wlanhb
     public :: la_wlanhe
     public :: la_wlanhf
     public :: la_wlanhp
     public :: la_wlanhs
     public :: la_wlanht
     public :: la_wlansb
     public :: la_wlansp
     public :: la_wlansy
     public :: la_wlantb
     public :: la_wlantp
     public :: la_wlantr

     contains

     !> SLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(sp) function la_slangb(norm,n,kl,ku,ab,ldab,work)
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
           la_slangb = value
           return
     end function la_slangb
     !> DLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(dp) function la_dlangb(norm,n,kl,ku,ab,ldab,work)
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
           la_dlangb = value
           return
     end function la_dlangb
     !> QLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(qp) function la_qlangb(norm,n,kl,ku,ab,ldab,work)
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
           la_qlangb = value
           return
     end function la_qlangb

     !> SLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(sp) function la_slange(norm,m,n,a,lda,work)
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
           la_slange = value
           return
     end function la_slange
     !> DLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(dp) function la_dlange(norm,m,n,a,lda,work)
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
           la_dlange = value
           return
     end function la_dlange
     !> QLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real matrix A.

     real(qp) function la_qlange(norm,m,n,a,lda,work)
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
           la_qlange = value
           return
     end function la_qlange

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

     !> SLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(sp) function la_slanhs(norm,n,a,lda,work)
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
           la_slanhs = value
           return
     end function la_slanhs
     !> DLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(dp) function la_dlanhs(norm,n,a,lda,work)
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
           la_dlanhs = value
           return
     end function la_dlanhs
     !> QLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(qp) function la_qlanhs(norm,n,a,lda,work)
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
           la_qlanhs = value
           return
     end function la_qlanhs

     !> SLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(sp) function la_slansb(norm,uplo,n,k,ab,ldab,work)
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
           la_slansb = value
           return
     end function la_slansb
     !> DLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(dp) function la_dlansb(norm,uplo,n,k,ab,ldab,work)
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
           la_dlansb = value
           return
     end function la_dlansb
     !> QLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(qp) function la_qlansb(norm,uplo,n,k,ab,ldab,work)
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
           la_qlansb = value
           return
     end function la_qlansb

     !> SLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(sp) function la_slansf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(sp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           if (n == 0) then
              la_slansf = zero
              return
           else if (n == 1) then
              la_slansf = abs(a(0))
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
           la_slansf = value
           return
     end function la_slansf
     !> DLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(dp) function la_dlansf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(dp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n == 0) then
              la_dlansf = zero
              return
           else if (n == 1) then
              la_dlansf = abs(a(0))
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
           la_dlansf = value
           return
     end function la_dlansf
     !> QLANSF: returns the value of the one norm, or the Frobenius norm, or
     !> the infinity norm, or the element of largest absolute value of a
     !> real symmetric matrix A in RFP format.

     real(qp) function la_qlansf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(qp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           if (n == 0) then
              la_qlansf = zero
              return
           else if (n == 1) then
              la_qlansf = abs(a(0))
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
           la_qlansf = value
           return
     end function la_qlansf

     !> SLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(sp) function la_slansp(norm,uplo,n,ap,work)
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
           la_slansp = value
           return
     end function la_slansp
     !> DLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(dp) function la_dlansp(norm,uplo,n,ap,work)
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
           la_dlansp = value
           return
     end function la_dlansp
     !> QLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A,  supplied in packed form.

     real(qp) function la_qlansp(norm,uplo,n,ap,work)
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
           la_qlansp = value
           return
     end function la_qlansp

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

     !> SLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(sp) function la_slansy(norm,uplo,n,a,lda,work)
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
           la_slansy = value
           return
     end function la_slansy
     !> DLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(dp) function la_dlansy(norm,uplo,n,a,lda,work)
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
           la_dlansy = value
           return
     end function la_dlansy
     !> QLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> real symmetric matrix A.

     real(qp) function la_qlansy(norm,uplo,n,a,lda,work)
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
           la_qlansy = value
           return
     end function la_qlansy

     !> SLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(sp) function la_slantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_slantb = value
           return
     end function la_slantb
     !> DLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(dp) function la_dlantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_dlantb = value
           return
     end function la_dlantb
     !> QLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(qp) function la_qlantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_qlantb = value
           return
     end function la_qlantb

     !> SLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(sp) function la_slantp(norm,uplo,diag,n,ap,work)
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
           la_slantp = value
           return
     end function la_slantp
     !> DLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(dp) function la_dlantp(norm,uplo,diag,n,ap,work)
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
           la_dlantp = value
           return
     end function la_dlantp
     !> QLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(qp) function la_qlantp(norm,uplo,diag,n,ap,work)
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
           la_qlantp = value
           return
     end function la_qlantp

     !> SLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(sp) function la_slantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_slantr = value
           return
     end function la_slantr
     !> DLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(dp) function la_dlantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_dlantr = value
           return
     end function la_dlantr
     !> QLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(qp) function la_qlantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_qlantr = value
           return
     end function la_qlantr

     !> CLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(sp) function la_clangb(norm,n,kl,ku,ab,ldab,work)
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
           la_clangb = value
           return
     end function la_clangb
     !> ZLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(dp) function la_zlangb(norm,n,kl,ku,ab,ldab,work)
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
           la_zlangb = value
           return
     end function la_zlangb
     !> WLANGB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n band matrix  A,  with kl sub-diagonals and ku super-diagonals.

     real(qp) function la_wlangb(norm,n,kl,ku,ab,ldab,work)
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
           la_wlangb = value
           return
     end function la_wlangb

     !> CLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(sp) function la_clange(norm,m,n,a,lda,work)
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
           la_clange = value
           return
     end function la_clange
     !> ZLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(dp) function la_zlange(norm,m,n,a,lda,work)
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
           la_zlange = value
           return
     end function la_zlange
     !> WLANGE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex matrix A.

     real(qp) function la_wlange(norm,m,n,a,lda,work)
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
           la_wlange = value
           return
     end function la_wlange

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

     !> CLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(sp) function la_clanhb(norm,uplo,n,k,ab,ldab,work)
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
           la_clanhb = value
           return
     end function la_clanhb
     !> ZLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(dp) function la_zlanhb(norm,uplo,n,k,ab,ldab,work)
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
           la_zlanhb = value
           return
     end function la_zlanhb
     !> WLANHB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n hermitian band matrix A,  with k super-diagonals.

     real(qp) function la_wlanhb(norm,uplo,n,k,ab,ldab,work)
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
           la_wlanhb = value
           return
     end function la_wlanhb

     !> CLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(sp) function la_clanhe(norm,uplo,n,a,lda,work)
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
           la_clanhe = value
           return
     end function la_clanhe
     !> ZLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(dp) function la_zlanhe(norm,uplo,n,a,lda,work)
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
           la_zlanhe = value
           return
     end function la_zlanhe
     !> WLANHE:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A.

     real(qp) function la_wlanhe(norm,uplo,n,a,lda,work)
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
           la_wlanhe = value
           return
     end function la_wlanhe

     !> CLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(sp) function la_clanhf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(sp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              la_clanhf = zero
              return
           else if (n == 1) then
              la_clanhf = abs(real(a(0),KIND=sp))
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
           la_clanhf = value
           return
     end function la_clanhf
     !> ZLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(dp) function la_zlanhf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(dp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              la_zlanhf = zero
              return
           else if (n == 1) then
              la_zlanhf = abs(real(a(0),KIND=dp))
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
           la_zlanhf = value
           return
     end function la_zlanhf
     !> WLANHF:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex Hermitian matrix A in RFP format.

     real(qp) function la_wlanhf(norm,transr,uplo,n,a,work)
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
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,ifm,ilu,noe,n1,k,l,lda
           real(qp) :: scale,s,value,aa,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,sqrt
           ! Executable Statements
           if (n == 0) then
              la_wlanhf = zero
              return
           else if (n == 1) then
              la_wlanhf = abs(real(a(0),KIND=qp))
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
           la_wlanhf = value
           return
     end function la_wlanhf

     !> CLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(sp) function la_clanhp(norm,uplo,n,ap,work)
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
           la_clanhp = value
           return
     end function la_clanhp
     !> ZLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(dp) function la_zlanhp(norm,uplo,n,ap,work)
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
           la_zlanhp = value
           return
     end function la_zlanhp
     !> WLANHP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex hermitian matrix A,  supplied in packed form.

     real(qp) function la_wlanhp(norm,uplo,n,ap,work)
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
           la_wlanhp = value
           return
     end function la_wlanhp

     !> CLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(sp) function la_clanhs(norm,n,a,lda,work)
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
           la_clanhs = value
           return
     end function la_clanhs
     !> ZLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(dp) function la_zlanhs(norm,n,a,lda,work)
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
           la_zlanhs = value
           return
     end function la_zlanhs
     !> WLANHS:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> Hessenberg matrix A.

     real(qp) function la_wlanhs(norm,n,a,lda,work)
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
           la_wlanhs = value
           return
     end function la_wlanhs

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

     !> CLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(sp) function la_clansb(norm,uplo,n,k,ab,ldab,work)
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
           la_clansb = value
           return
     end function la_clansb
     !> ZLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(dp) function la_zlansb(norm,uplo,n,k,ab,ldab,work)
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
           la_zlansb = value
           return
     end function la_zlansb
     !> WLANSB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n symmetric band matrix A,  with k super-diagonals.

     real(qp) function la_wlansb(norm,uplo,n,k,ab,ldab,work)
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
           la_wlansb = value
           return
     end function la_wlansb

     !> CLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(sp) function la_clansp(norm,uplo,n,ap,work)
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
           la_clansp = value
           return
     end function la_clansp
     !> ZLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(dp) function la_zlansp(norm,uplo,n,ap,work)
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
           la_zlansp = value
           return
     end function la_zlansp
     !> WLANSP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A,  supplied in packed form.

     real(qp) function la_wlansp(norm,uplo,n,ap,work)
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
           la_wlansp = value
           return
     end function la_wlansp

     !> CLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(sp) function la_clansy(norm,uplo,n,a,lda,work)
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
           la_clansy = value
           return
     end function la_clansy
     !> ZLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(dp) function la_zlansy(norm,uplo,n,a,lda,work)
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
           la_zlansy = value
           return
     end function la_zlansy
     !> WLANSY:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> complex symmetric matrix A.

     real(qp) function la_wlansy(norm,uplo,n,a,lda,work)
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
           la_wlansy = value
           return
     end function la_wlansy

     !> CLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(sp) function la_clantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_clantb = value
           return
     end function la_clantb
     !> ZLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(dp) function la_zlantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_zlantb = value
           return
     end function la_zlantb
     !> WLANTB:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the element of  largest absolute value  of an
     !> n by n triangular band matrix A,  with ( k + 1 ) diagonals.

     real(qp) function la_wlantb(norm,uplo,diag,n,k,ab,ldab,work)
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
           la_wlantb = value
           return
     end function la_wlantb

     !> CLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(sp) function la_clantp(norm,uplo,diag,n,ap,work)
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
           la_clantp = value
           return
     end function la_clantp
     !> ZLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(dp) function la_zlantp(norm,uplo,diag,n,ap,work)
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
           la_zlantp = value
           return
     end function la_zlantp
     !> WLANTP:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> triangular matrix A, supplied in packed form.

     real(qp) function la_wlantp(norm,uplo,diag,n,ap,work)
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
           la_wlantp = value
           return
     end function la_wlantp

     !> CLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(sp) function la_clantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_clantr = value
           return
     end function la_clantr
     !> ZLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(dp) function la_zlantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_zlantr = value
           return
     end function la_zlantr
     !> WLANTR:  returns the value of the one norm,  or the Frobenius norm, or
     !> the  infinity norm,  or the  element of  largest absolute value  of a
     !> trapezoidal or triangular matrix A.

     real(qp) function la_wlantr(norm,uplo,diag,m,n,a,lda,work)
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
           la_wlantr = value
           return
     end function la_wlantr

end module la_lapack_blas_like_mnorm
