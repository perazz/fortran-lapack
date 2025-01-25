module la_linalg_blas_aux
     use la_linalg_constants
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_dcabs1
     public :: la_icamax
     public :: la_idamax
     public :: la_isamax
     public :: la_izamax
     public :: la_lsame
     public :: la_scabs1
     public :: la_xerbla
     public :: la_xerbla_array
     public :: la_qcabs1
     public :: la_iqamax
     public :: la_iwamax

     contains

     !> DCABS1: computes |Re(.)| + |Im(.)| of a double complex number

     pure real(dp) function la_dcabs1(z)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: z
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag
           la_dcabs1 = abs(real(z,KIND=dp)) + abs(aimag(z))
           return
     end function la_dcabs1

     !> ISAMAX: finds the index of the first element having maximum absolute value.

     pure integer(ilp) function la_isamax(n,sx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: smax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           la_isamax = 0
           if (n < 1 .or. incx <= 0) return
           la_isamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              smax = abs(sx(1))
              do i = 2,n
                 if (abs(sx(i)) > smax) then
                    la_isamax = i
                    smax = abs(sx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              smax = abs(sx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(sx(ix)) > smax) then
                    la_isamax = i
                    smax = abs(sx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_isamax

     !> IZAMAX: finds the index of the first element having maximum |Re(.)| + |Im(.)|

     pure integer(ilp) function la_izamax(n,zx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dmax
           integer(ilp) :: i,ix
           la_izamax = 0
           if (n < 1 .or. incx <= 0) return
           la_izamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = la_dcabs1(zx(1))
              do i = 2,n
                 if (la_dcabs1(zx(i)) > dmax) then
                    la_izamax = i
                    dmax = la_dcabs1(zx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = la_dcabs1(zx(1))
              ix = ix + incx
              do i = 2,n
                 if (la_dcabs1(zx(ix)) > dmax) then
                    la_izamax = i
                    dmax = la_dcabs1(zx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_izamax

     !> LSAME: returns .TRUE. if CA is the same letter as CB regardless of
     !> case.

     pure logical(lk) function la_lsame(ca,cb)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: ca,cb
       ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: ichar
           ! Local Scalars
           integer(ilp) :: inta,intb,zcode
           ! test if the characters are equal
           la_lsame = ca == cb
           if (la_lsame) return
           ! now test for equivalence if both characters are alphabetic.
           zcode = ichar('Z')
           ! use 'z' rather than 'a' so that ascii can be detected on prime
           ! machines, on which ichar returns a value with bit 8 set.
           ! ichar('a') on prime machines returns 193 which is the same as
           ! ichar('a') on an ebcdic machine.
           inta = ichar(ca)
           intb = ichar(cb)
           if (zcode == 90 .or. zcode == 122) then
              ! ascii is assumed - zcode is the ascii code of either lower or
              ! upper case 'z'.
               if (inta >= 97 .and. inta <= 122) inta = inta - 32
               if (intb >= 97 .and. intb <= 122) intb = intb - 32
           else if (zcode == 233 .or. zcode == 169) then
              ! ebcdic is assumed - zcode is the ebcdic code of either lower or
              ! upper case 'z'.
               if (inta >= 129 .and. inta <= 137 .or. inta >= 145 .and. inta <= 153 .or. inta >= 162 .and. &
                         inta <= 169) inta = inta + 64
               if (intb >= 129 .and. intb <= 137 .or. intb >= 145 .and. intb <= 153 .or. intb >= 162 .and. &
                         intb <= 169) intb = intb + 64
           else if (zcode == 218 .or. zcode == 250) then
              ! ascii is assumed, on prime machines - zcode is the ascii code
              ! plus 128 of either lower or upper case 'z'.
               if (inta >= 225 .and. inta <= 250) inta = inta - 32
               if (intb >= 225 .and. intb <= 250) intb = intb - 32
           end if
           la_lsame = inta == intb
           ! return
     end function la_lsame

     !> SCABS1: computes |Re(.)| + |Im(.)| of a complex number

     pure real(sp) function la_scabs1(z)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: z
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: abs,aimag,real
           la_scabs1 = abs(real(z,KIND=sp)) + abs(aimag(z))
           return
     end function la_scabs1

     !> XERBLA:  is an error handler for the LAPACK routines.
     !> It is called by an LAPACK routine if an input parameter has an
     !> invalid value.  A message is printed and execution stops.
     !> Installers may consider modifying the STOP statement in order to
     !> call system-specific exception-handling facilities.

     pure subroutine la_xerbla(srname,info)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character(len=*),intent(in) :: srname
           integer(ilp),intent(in) :: info
       ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: len_trim
           ! Executable Statements
9999  format(' ** ON ENTRY TO ',a,' PARAMETER NUMBER ',i2,' HAD ','AN ILLEGAL VALUE')
                
     end subroutine la_xerbla

     !> XERBLA_ARRAY: assists other languages in calling XERBLA, the LAPACK
     !> and BLAS error handler.  Rather than taking a Fortran string argument
     !> as the function's name, XERBLA_ARRAY takes an array of single
     !> characters along with the array's length.  XERBLA_ARRAY then copies
     !> up to 32 characters of that array into a Fortran string and passes
     !> that to XERBLA.  If called with a non-positive SRNAME_LEN,
     !> XERBLA_ARRAY will call XERBLA with a string of all blank characters.
     !> Say some macro or other device makes XERBLA_ARRAY available to C99
     !> by a name lapack_xerbla and with a common Fortran calling convention.
     !> Then a C99 program could invoke XERBLA via:
     !> {
     !> int flen = strlen(__func__);
     !> lapack_xerbla(__func__,
     !> }
     !> Providing XERBLA_ARRAY is not necessary for intercepting LAPACK
     !> errors.  XERBLA_ARRAY calls XERBLA.

     pure subroutine la_xerbla_array(srname_array,srname_len,info)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: srname_len,info
           ! Array Arguments
           character(1),intent(in) :: srname_array(srname_len)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i
           ! Local Arrays
           character*32 srname
           ! Intrinsic Functions
           intrinsic :: min,len
           ! Executable Statements
           srname = ''
           do i = 1,min(srname_len,len(srname))
              srname(i:i) = srname_array(i)
           end do
           call la_xerbla(srname,info)
           return
     end subroutine la_xerbla_array

     !> DCABS1: computes |Re(.)| + |Im(.)| of a double complex number

     pure real(qp) function la_qcabs1(z)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: z
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag
           la_qcabs1 = abs(real(z,KIND=qp)) + abs(aimag(z))
           return
     end function la_qcabs1

     !> IDAMAX: finds the index of the first element having maximum absolute value.

     pure integer(ilp) function la_iqamax(n,dx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(qp),intent(in) :: dx(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: dmax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           la_iqamax = 0
           if (n < 1 .or. incx <= 0) return
           la_iqamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = abs(dx(1))
              do i = 2,n
                 if (abs(dx(i)) > dmax) then
                    la_iqamax = i
                    dmax = abs(dx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = abs(dx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(dx(ix)) > dmax) then
                    la_iqamax = i
                    dmax = abs(dx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_iqamax

     !> IZAMAX: finds the index of the first element having maximum |Re(.)| + |Im(.)|

     pure integer(ilp) function la_iwamax(n,zx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(in) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: dmax
           integer(ilp) :: i,ix
           la_iwamax = 0
           if (n < 1 .or. incx <= 0) return
           la_iwamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = la_qcabs1(zx(1))
              do i = 2,n
                 if (la_qcabs1(zx(i)) > dmax) then
                    la_iwamax = i
                    dmax = la_qcabs1(zx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = la_qcabs1(zx(1))
              ix = ix + incx
              do i = 2,n
                 if (la_qcabs1(zx(ix)) > dmax) then
                    la_iwamax = i
                    dmax = la_qcabs1(zx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_iwamax

     !> ICAMAX: finds the index of the first element having maximum |Re(.)| + |Im(.)|

     pure integer(ilp) function la_icamax(n,cx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: smax
           integer(ilp) :: i,ix
           la_icamax = 0
           if (n < 1 .or. incx <= 0) return
           la_icamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              smax = la_scabs1(cx(1))
              do i = 2,n
                 if (la_scabs1(cx(i)) > smax) then
                    la_icamax = i
                    smax = la_scabs1(cx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              smax = la_scabs1(cx(1))
              ix = ix + incx
              do i = 2,n
                 if (la_scabs1(cx(ix)) > smax) then
                    la_icamax = i
                    smax = la_scabs1(cx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_icamax

     !> IDAMAX: finds the index of the first element having maximum absolute value.

     pure integer(ilp) function la_idamax(n,dx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(dp),intent(in) :: dx(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dmax
           integer(ilp) :: i,ix
           ! Intrinsic Functions
           intrinsic :: abs
           la_idamax = 0
           if (n < 1 .or. incx <= 0) return
           la_idamax = 1
           if (n == 1) return
           if (incx == 1) then
              ! code for increment equal to 1
              dmax = abs(dx(1))
              do i = 2,n
                 if (abs(dx(i)) > dmax) then
                    la_idamax = i
                    dmax = abs(dx(i))
                 end if
              end do
           else
              ! code for increment not equal to 1
              ix = 1
              dmax = abs(dx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(dx(ix)) > dmax) then
                    la_idamax = i
                    dmax = abs(dx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     end function la_idamax

end module la_linalg_blas_aux
