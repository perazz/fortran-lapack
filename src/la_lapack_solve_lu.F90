!> LU drivers: general, banded and tridiagonal systems
module la_lapack_solve_lu
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level3_gen
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_mnorm
     use la_lapack_solve_lu_comp
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgtsv
     public :: la_sgbsv
     public :: la_sgbsvx
     public :: la_sgtsvx
     public :: la_sgesv
     public :: la_sgesvx
     public :: la_dgtsv
     public :: la_dgbsv
     public :: la_dgbsvx
     public :: la_dgtsvx
     public :: la_dsgesv
     public :: la_dgesv
     public :: la_dgesvx
#ifdef LA_WITH_XDP
     public :: la_xgtsv
     public :: la_xgbsv
     public :: la_xgbsvx
     public :: la_xgtsvx
     public :: la_xdgesv
     public :: la_xgesv
     public :: la_xgesvx
#endif
#ifdef LA_WITH_QP
     public :: la_qgtsv
     public :: la_qgbsv
     public :: la_qgbsvx
     public :: la_qgtsvx
     public :: la_qdgesv
     public :: la_qgesv
     public :: la_qgesvx
#endif
     public :: la_cgtsv
     public :: la_cgbsv
     public :: la_cgbsvx
     public :: la_cgtsvx
     public :: la_cgesv
     public :: la_cgesvx
     public :: la_zgtsv
     public :: la_zgbsv
     public :: la_zgbsvx
     public :: la_zgtsvx
     public :: la_zcgesv
     public :: la_zgesv
     public :: la_zgesvx
#ifdef LA_WITH_XDP
     public :: la_ygtsv
     public :: la_ygbsv
     public :: la_ygbsvx
     public :: la_ygtsvx
     public :: la_yzgesv
     public :: la_ygesv
     public :: la_ygesvx
#endif
#ifdef LA_WITH_QP
     public :: la_wgtsv
     public :: la_wgbsv
     public :: la_wgbsvx
     public :: la_wgtsvx
     public :: la_wzgesv
     public :: la_wgesv
     public :: la_wgesvx
#endif

     contains

     !> SGTSV:  solves the equation
     !> A*X = B,
     !> where A is an n by n tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T*X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_sgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: fact,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SGTSV ',-info)
              return
           end if
           if (n == 0) return
           if (nrhs == 1) then
              loop_10: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end do loop_10
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           else
              loop_40: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end do loop_40
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           end if
           ! back solve with the matrix u from the factorization.
           if (nrhs <= 2) then
              j = 1
              70 continue
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do i = n - 2,1,-1
                 b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

              end do
              if (j < nrhs) then
                 j = j + 1
                 go to 70
              end if
           else
              do j = 1,nrhs
                 b(n,j) = b(n,j)/d(n)
                 if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
                 do i = n - 2,1,-1
                    b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

                 end do
              end do
           end if
           return
     end subroutine la_sgtsv
     !> DGTSV:  solves the equation
     !> A*X = B,
     !> where A is an n by n tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T*X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_dgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: fact,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DGTSV ',-info)
              return
           end if
           if (n == 0) return
           if (nrhs == 1) then
              loop_10: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end do loop_10
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           else
              loop_40: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end do loop_40
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           end if
           ! back solve with the matrix u from the factorization.
           if (nrhs <= 2) then
              j = 1
              70 continue
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do i = n - 2,1,-1
                 b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

              end do
              if (j < nrhs) then
                 j = j + 1
                 go to 70
              end if
           else
              do j = 1,nrhs
                 b(n,j) = b(n,j)/d(n)
                 if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
                 do i = n - 2,1,-1
                    b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

                 end do
              end do
           end if
           return
     end subroutine la_dgtsv
#ifdef LA_WITH_XDP
     !> XGTSV:  solves the equation
     !> A*X = B,
     !> where A is an n by n tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T*X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_xgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(xdp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: fact,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XGTSV ',-info)
              return
           end if
           if (n == 0) return
           if (nrhs == 1) then
              loop_10: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end do loop_10
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           else
              loop_40: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end do loop_40
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           end if
           ! back solve with the matrix u from the factorization.
           if (nrhs <= 2) then
              j = 1
              70 continue
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do i = n - 2,1,-1
                 b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

              end do
              if (j < nrhs) then
                 j = j + 1
                 go to 70
              end if
           else
              do j = 1,nrhs
                 b(n,j) = b(n,j)/d(n)
                 if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
                 do i = n - 2,1,-1
                    b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

                 end do
              end do
           end if
           return
     end subroutine la_xgtsv
#endif
#ifdef LA_WITH_QP
     !> QGTSV:  solves the equation
     !> A*X = B,
     !> where A is an n by n tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T*X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_qgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: fact,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QGTSV ',-info)
              return
           end if
           if (n == 0) return
           if (nrhs == 1) then
              loop_10: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end do loop_10
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       b(i + 1,1) = b(i + 1,1) - fact*b(i,1)
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    temp = b(i,1)
                    b(i,1) = b(i + 1,1)
                    b(i + 1,1) = temp - fact*b(i + 1,1)
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           else
              loop_40: do i = 1,n - 2
                 if (abs(d(i)) >= abs(dl(i))) then
                    ! no row interchange required
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                    dl(i) = zero
                 else
                    ! interchange rows i and i+1
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    dl(i) = du(i + 1)
                    du(i + 1) = -fact*dl(i)
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end do loop_40
              if (n > 1) then
                 i = n - 1
                 if (abs(d(i)) >= abs(dl(i))) then
                    if (d(i) /= zero) then
                       fact = dl(i)/d(i)
                       d(i + 1) = d(i + 1) - fact*du(i)
                       do j = 1,nrhs
                          b(i + 1,j) = b(i + 1,j) - fact*b(i,j)
                       end do
                    else
                       info = i
                       return
                    end if
                 else
                    fact = d(i)/dl(i)
                    d(i) = dl(i)
                    temp = d(i + 1)
                    d(i + 1) = du(i) - fact*temp
                    du(i) = temp
                    do j = 1,nrhs
                       temp = b(i,j)
                       b(i,j) = b(i + 1,j)
                       b(i + 1,j) = temp - fact*b(i + 1,j)
                    end do
                 end if
              end if
              if (d(n) == zero) then
                 info = n
                 return
              end if
           end if
           ! back solve with the matrix u from the factorization.
           if (nrhs <= 2) then
              j = 1
              70 continue
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do i = n - 2,1,-1
                 b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

              end do
              if (j < nrhs) then
                 j = j + 1
                 go to 70
              end if
           else
              do j = 1,nrhs
                 b(n,j) = b(n,j)/d(n)
                 if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
                 do i = n - 2,1,-1
                    b(i,j) = (b(i,j) - du(i)*b(i + 1,j) - dl(i)*b(i + 2,j))/d(i)

                 end do
              end do
           end if
           return
     end subroutine la_qgtsv
#endif

     !> SGBSV: computes the solution to a real system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_sgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_sgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_sgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_sgbsv
     !> DGBSV: computes the solution to a real system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_dgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_dgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_dgbsv
#ifdef LA_WITH_XDP
     !> XGBSV: computes the solution to a real system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_xgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_xgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_xgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_xgbsv
#endif
#ifdef LA_WITH_QP
     !> QGBSV: computes the solution to a real system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_qgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_qgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_qgbsv
#endif

     !> SGBSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_sgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),c(*),r(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================
        ! moved setting of info = n+1 so info does not subsequently get
        ! overwritten.  sven, 17 mar 05.
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(sp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_sgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_slaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_scopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_sgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_slantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_slangb(norm,n,kl,ku,ab,ldab,work)
           rpvgrw = la_slantb('M','U','N',n,kl + ku,afb,ldafb,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_slangb('M',n,kl,ku,ab,ldab,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_sgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution matrix x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_sgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_sgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           work(1) = rpvgrw
           return
     end subroutine la_sgbsvx
     !> DGBSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),c(*),r(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(dp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_dgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_dlaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_dcopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_dgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_dlantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_dlangb(norm,n,kl,ku,ab,ldab,work)
           rpvgrw = la_dlantb('M','U','N',n,kl + ku,afb,ldafb,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_dlangb('M',n,kl,ku,ab,ldab,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_dgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution matrix x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_dgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           work(1) = rpvgrw
           return
     end subroutine la_dgbsvx
#ifdef LA_WITH_XDP
     !> XGBSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_xgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),c(*),r(*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(xdp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_xlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_xgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_xlaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_xcopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_xgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_xlantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_xlangb(norm,n,kl,ku,ab,ldab,work)
           rpvgrw = la_xlantb('M','U','N',n,kl + ku,afb,ldafb,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_xlangb('M',n,kl,ku,ab,ldab,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_xgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution matrix x.
           call la_xlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_xgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_xgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           work(1) = rpvgrw
           return
     end subroutine la_xgbsvx
#endif
#ifdef LA_WITH_QP
     !> QGBSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),c(*),r(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(qp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_qgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_qlaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_qcopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_qgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_qlantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_qlangb(norm,n,kl,ku,ab,ldab,work)
           rpvgrw = la_qlantb('M','U','N',n,kl + ku,afb,ldafb,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_qlangb('M',n,kl,ku,ab,ldab,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_qgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution matrix x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_qgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           work(1) = rpvgrw
           return
     end subroutine la_qgbsvx
#endif

     !> SGTSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B or A**T * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_sgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(sp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('SGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_scopy(n,d,1,df,1)
              if (n > 1) then
                 call la_scopy(n - 1,dl,1,dlf,1)
                 call la_scopy(n - 1,du,1,duf,1)
              end if
              call la_sgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_slangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_sgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution vectors x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_sgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_sgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_sgtsvx
     !> DGTSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B or A**T * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_dgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(dp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_dcopy(n,d,1,df,1)
              if (n > 1) then
                 call la_dcopy(n - 1,dl,1,dlf,1)
                 call la_dcopy(n - 1,du,1,duf,1)
              end if
              call la_dgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_dlangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_dgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution vectors x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_dgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dgtsvx
#ifdef LA_WITH_XDP
     !> XGTSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B or A**T * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_xgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(xdp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('XGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_xcopy(n,d,1,df,1)
              if (n > 1) then
                 call la_xcopy(n - 1,dl,1,dlf,1)
                 call la_xcopy(n - 1,du,1,duf,1)
              end if
              call la_xgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_xlangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_xgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution vectors x.
           call la_xlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_xgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_xgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_xgtsvx
#endif
#ifdef LA_WITH_QP
     !> QGTSVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations A * X = B or A**T * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_qgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(qp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_qcopy(n,d,1,df,1)
              if (n > 1) then
                 call la_qcopy(n - 1,dl,1,dlf,1)
                 call la_qcopy(n - 1,du,1,duf,1)
              end if
              call la_qgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_qlangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_qgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,iwork,info)

           ! compute the solution vectors x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_qgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qgtsvx
#endif

     !> DSGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> DSGESV first attempts to factorize the matrix in SINGLE PRECISION
     !> and use this factorization within an iterative refinement procedure
     !> to produce a solution with DOUBLE PRECISION normwise backward error
     !> quality (see below). If the approach fails the method switches to a
     !> DOUBLE PRECISION factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio SINGLE PRECISION performance over DOUBLE PRECISION
     !> performance is too small. A reasonable strategy should take the
     !> number of right-hand sides and the size of the matrix into account.
     !> This might be done with a call to ILAENV in the future. Up to now, we
     !> always try iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by DLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_dsgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,iter,info)
        use la_constants_dp,only:negone,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(out) :: swork(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: b(ldb,*)
           real(dp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(dp),parameter :: bwdmax = 1.0e+00_dp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(dp) :: anrm,cte,eps,rnrm,xnrm
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DSGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip single precision iterative refinement if a priori slower
           ! than double precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_dlange('I',n,n,a,lda,work)
           eps = la_dlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=dp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from double precision to single precision and store the
           ! result in sx.
           call la_dlag2s(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from double precision to single precision and store the
           ! result in sa.
           call la_dlag2s(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_sgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_sgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to double precision
           call la_slag2d(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_dlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x,ldx, &
                     one,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = abs(x(la_idamax(n,x(1,i),1),i))
              rnrm = abs(work(la_idamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from double precision to single precision
              ! and store the result in sx.
              call la_dlag2s(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_sgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to double precision and update the current
              ! iterate.
              call la_slag2d(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_daxpy(n,one,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_dlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x, &
                        ldx,one,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = abs(x(la_idamax(n,x(1,i),1),i))
                 rnrm = abs(work(la_idamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow up
           ! on double precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to double precision.
           call la_dgetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_dlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_dgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_dsgesv
#ifdef LA_WITH_XDP
     !> XDGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> XDGESV first attempts to factorize the matrix in SINGLE PRECISION
     !> and use this factorization within an iterative refinement procedure
     !> to produce a solution with EXTENDED PRECISION normwise backward error
     !> quality (see below). If the approach fails the method switches to a
     !> EXTENDED PRECISION factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio SINGLE PRECISION performance over EXTENDED PRECISION
     !> performance is too small. A reasonable strategy should take the
     !> number of right-hand sides and the size of the matrix into account.
     !> This might be done with a call to ILAENV in the future. Up to now, we
     !> always try iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by XLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_xdgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,iter,info)
        use la_constants_xdp,only:negone,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(out) :: swork(*)
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: b(ldb,*)
           real(xdp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(xdp),parameter :: bwdmax = 1.0e+00_xdp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(xdp) :: anrm,cte,eps,rnrm,xnrm
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XDGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than extended precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_xlange('I',n,n,a,lda,work)
           eps = la_xlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=xdp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from extended precision to double precision and store the
           ! result in sx.
           call la_xlag2d(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from extended precision to double precision and store the
           ! result in sa.
           call la_xlag2d(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_dgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_dgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to extended precision
           call la_dlag2x(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_xlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x,ldx, &
                     one,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = abs(x(la_ixamax(n,x(1,i),1),i))
              rnrm = abs(work(la_ixamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from extended precision to double precision
              ! and store the result in sx.
              call la_xlag2d(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_dgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to extended precision and update the current
              ! iterate.
              call la_dlag2x(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_xaxpy(n,one,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_xlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x, &
                        ldx,one,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = abs(x(la_ixamax(n,x(1,i),1),i))
                 rnrm = abs(work(la_ixamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow up
           ! on extended precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to extended precision.
           call la_xgetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_xlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_xgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_xdgesv
#endif
#ifdef LA_WITH_QP
     !> QDGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> QDGESV first attempts to factorize the matrix in SINGLE PRECISION
     !> and use this factorization within an iterative refinement procedure
     !> to produce a solution with QUAD PRECISION normwise backward error
     !> quality (see below). If the approach fails the method switches to a
     !> QUAD PRECISION factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio SINGLE PRECISION performance over QUAD PRECISION
     !> performance is too small. A reasonable strategy should take the
     !> number of right-hand sides and the size of the matrix into account.
     !> This might be done with a call to ILAENV in the future. Up to now, we
     !> always try iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by QLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_qdgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,iter,info)
        use la_constants_qp,only:negone,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(out) :: swork(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: b(ldb,*)
           real(qp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(qp),parameter :: bwdmax = 1.0e+00_qp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(qp) :: anrm,cte,eps,rnrm,xnrm
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QDGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than quad precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_qlange('I',n,n,a,lda,work)
           eps = la_qlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=qp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from quad precision to double precision and store the
           ! result in sx.
           call la_qlag2d(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from quad precision to double precision and store the
           ! result in sa.
           call la_qlag2d(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_dgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_dgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to quad precision
           call la_dlag2q(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_qlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x,ldx, &
                     one,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = abs(x(la_iqamax(n,x(1,i),1),i))
              rnrm = abs(work(la_iqamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from quad precision to double precision
              ! and store the result in sx.
              call la_qlag2d(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_dgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to quad precision and update the current
              ! iterate.
              call la_dlag2q(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_qaxpy(n,one,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_qlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,negone,a,lda,x, &
                        ldx,one,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = abs(x(la_iqamax(n,x(1,i),1),i))
                 rnrm = abs(work(la_iqamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow up
           ! on quad precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to quad precision.
           call la_qgetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_qlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_qgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_qdgesv
#endif

     !> SGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_sgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_sgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_sgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_sgesv
     !> DGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_dgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_dgesv
#ifdef LA_WITH_XDP
     !> XGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_xgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_xgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_xgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_xgesv
#endif
#ifdef LA_WITH_QP
     !> QGESV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_qgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_qgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_qgesv
#endif

     !> SGESVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_sgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),c(*),r(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_sgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_slaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_slacpy('FULL',n,n,a,lda,af,ldaf)
              call la_sgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_slantr('M','U','N',info,info,af,ldaf,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_slange('M',n,info,a,lda,work)/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_slange(norm,n,n,a,lda,work)
           rpvgrw = la_slantr('M','U','N',n,n,af,ldaf,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_slange('M',n,n,a,lda,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_sgecon(norm,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_sgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_sgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           work(1) = rpvgrw
           return
     end subroutine la_sgesvx
     !> DGESVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),c(*),r(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_dgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_dlaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_dlacpy('FULL',n,n,a,lda,af,ldaf)
              call la_dgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_dlantr('M','U','N',info,info,af,ldaf,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_dlange('M',n,info,a,lda,work)/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_dlange(norm,n,n,a,lda,work)
           rpvgrw = la_dlantr('M','U','N',n,n,af,ldaf,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_dlange('M',n,n,a,lda,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_dgecon(norm,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_dgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           work(1) = rpvgrw
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dgesvx
#ifdef LA_WITH_XDP
     !> XGESVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_xgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),c(*),r(*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(xdp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_xlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_xgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_xlaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_xlacpy('FULL',n,n,a,lda,af,ldaf)
              call la_xgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_xlantr('M','U','N',info,info,af,ldaf,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_xlange('M',n,info,a,lda,work)/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_xlange(norm,n,n,a,lda,work)
           rpvgrw = la_xlantr('M','U','N',n,n,af,ldaf,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_xlange('M',n,n,a,lda,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_xgecon(norm,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_xlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_xgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_xgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           work(1) = rpvgrw
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_xgesvx
#endif
#ifdef LA_WITH_QP
     !> QGESVX: uses the LU factorization to compute the solution to a real
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),c(*),r(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_qgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_qlaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_qlacpy('FULL',n,n,a,lda,af,ldaf)
              call la_qgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_qlantr('M','U','N',info,info,af,ldaf,work)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_qlange('M',n,info,a,lda,work)/rpvgrw
                 end if
                 work(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_qlange(norm,n,n,a,lda,work)
           rpvgrw = la_qlantr('M','U','N',n,n,af,ldaf,work)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_qlange('M',n,n,a,lda,work)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_qgecon(norm,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_qgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           work(1) = rpvgrw
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qgesvx
#endif

     !> CGTSV:  solves the equation
     !> A*X = B,
     !> where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T *X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_cgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k
           complex(sp) :: mult,temp,zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CGTSV ',-info)
              return
           end if
           if (n == 0) return
           loop_30: do k = 1,n - 1
              if (dl(k) == czero) then
                 ! subdiagonal is czero, no elimination is required.
                 if (d(k) == czero) then
                    ! diagonal is czero: set info = k and return; a unique
                    ! solution can not be found.
                    info = k
                    return
                 end if
              else if (cabs1(d(k)) >= cabs1(dl(k))) then
                 ! no row interchange required
                 mult = dl(k)/d(k)
                 d(k + 1) = d(k + 1) - mult*du(k)
                 do j = 1,nrhs
                    b(k + 1,j) = b(k + 1,j) - mult*b(k,j)
                 end do
                 if (k < (n - 1)) dl(k) = czero
              else
                 ! interchange rows k and k+1
                 mult = d(k)/dl(k)
                 d(k) = dl(k)
                 temp = d(k + 1)
                 d(k + 1) = du(k) - mult*temp
                 if (k < (n - 1)) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                 end if
                 du(k) = temp
                 do j = 1,nrhs
                    temp = b(k,j)
                    b(k,j) = b(k + 1,j)
                    b(k + 1,j) = temp - mult*b(k + 1,j)
                 end do
              end if
           end do loop_30
           if (d(n) == czero) then
              info = n
              return
           end if
           ! back solve with the matrix u from the factorization.
           do j = 1,nrhs
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do k = n - 2,1,-1
                 b(k,j) = (b(k,j) - du(k)*b(k + 1,j) - dl(k)*b(k + 2,j))/d(k)

              end do
           end do
           return
     end subroutine la_cgtsv
     !> ZGTSV:  solves the equation
     !> A*X = B,
     !> where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T *X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_zgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k
           complex(dp) :: mult,temp,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZGTSV ',-info)
              return
           end if
           if (n == 0) return
           loop_30: do k = 1,n - 1
              if (dl(k) == czero) then
                 ! subdiagonal is czero, no elimination is required.
                 if (d(k) == czero) then
                    ! diagonal is czero: set info = k and return; a unique
                    ! solution can not be found.
                    info = k
                    return
                 end if
              else if (cabs1(d(k)) >= cabs1(dl(k))) then
                 ! no row interchange required
                 mult = dl(k)/d(k)
                 d(k + 1) = d(k + 1) - mult*du(k)
                 do j = 1,nrhs
                    b(k + 1,j) = b(k + 1,j) - mult*b(k,j)
                 end do
                 if (k < (n - 1)) dl(k) = czero
              else
                 ! interchange rows k and k+1
                 mult = d(k)/dl(k)
                 d(k) = dl(k)
                 temp = d(k + 1)
                 d(k + 1) = du(k) - mult*temp
                 if (k < (n - 1)) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                 end if
                 du(k) = temp
                 do j = 1,nrhs
                    temp = b(k,j)
                    b(k,j) = b(k + 1,j)
                    b(k + 1,j) = temp - mult*b(k + 1,j)
                 end do
              end if
           end do loop_30
           if (d(n) == czero) then
              info = n
              return
           end if
           ! back solve with the matrix u from the factorization.
           do j = 1,nrhs
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do k = n - 2,1,-1
                 b(k,j) = (b(k,j) - du(k)*b(k + 1,j) - dl(k)*b(k + 2,j))/d(k)

              end do
           end do
           return
     end subroutine la_zgtsv
#ifdef LA_WITH_XDP
     !> YGTSV:  solves the equation
     !> A*X = B,
     !> where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T *X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_ygtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(xdp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k
           complex(xdp) :: mult,temp,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YGTSV ',-info)
              return
           end if
           if (n == 0) return
           loop_30: do k = 1,n - 1
              if (dl(k) == czero) then
                 ! subdiagonal is czero, no elimination is required.
                 if (d(k) == czero) then
                    ! diagonal is czero: set info = k and return; a unique
                    ! solution can not be found.
                    info = k
                    return
                 end if
              else if (cabs1(d(k)) >= cabs1(dl(k))) then
                 ! no row interchange required
                 mult = dl(k)/d(k)
                 d(k + 1) = d(k + 1) - mult*du(k)
                 do j = 1,nrhs
                    b(k + 1,j) = b(k + 1,j) - mult*b(k,j)
                 end do
                 if (k < (n - 1)) dl(k) = czero
              else
                 ! interchange rows k and k+1
                 mult = d(k)/dl(k)
                 d(k) = dl(k)
                 temp = d(k + 1)
                 d(k + 1) = du(k) - mult*temp
                 if (k < (n - 1)) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                 end if
                 du(k) = temp
                 do j = 1,nrhs
                    temp = b(k,j)
                    b(k,j) = b(k + 1,j)
                    b(k + 1,j) = temp - mult*b(k + 1,j)
                 end do
              end if
           end do loop_30
           if (d(n) == czero) then
              info = n
              return
           end if
           ! back solve with the matrix u from the factorization.
           do j = 1,nrhs
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do k = n - 2,1,-1
                 b(k,j) = (b(k,j) - du(k)*b(k + 1,j) - dl(k)*b(k + 2,j))/d(k)

              end do
           end do
           return
     end subroutine la_ygtsv
#endif
#ifdef LA_WITH_QP
     !> WGTSV:  solves the equation
     !> A*X = B,
     !> where A is an N-by-N tridiagonal matrix, by Gaussian elimination with
     !> partial pivoting.
     !> Note that the equation  A**T *X = B  may be solved by interchanging the
     !> order of the arguments DU and DL.

     pure subroutine la_wgtsv(n,nrhs,dl,d,du,b,ldb,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: b(ldb,*),d(*),dl(*),du(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k
           complex(qp) :: mult,temp,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WGTSV ',-info)
              return
           end if
           if (n == 0) return
           loop_30: do k = 1,n - 1
              if (dl(k) == czero) then
                 ! subdiagonal is czero, no elimination is required.
                 if (d(k) == czero) then
                    ! diagonal is czero: set info = k and return; a unique
                    ! solution can not be found.
                    info = k
                    return
                 end if
              else if (cabs1(d(k)) >= cabs1(dl(k))) then
                 ! no row interchange required
                 mult = dl(k)/d(k)
                 d(k + 1) = d(k + 1) - mult*du(k)
                 do j = 1,nrhs
                    b(k + 1,j) = b(k + 1,j) - mult*b(k,j)
                 end do
                 if (k < (n - 1)) dl(k) = czero
              else
                 ! interchange rows k and k+1
                 mult = d(k)/dl(k)
                 d(k) = dl(k)
                 temp = d(k + 1)
                 d(k + 1) = du(k) - mult*temp
                 if (k < (n - 1)) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                 end if
                 du(k) = temp
                 do j = 1,nrhs
                    temp = b(k,j)
                    b(k,j) = b(k + 1,j)
                    b(k + 1,j) = temp - mult*b(k + 1,j)
                 end do
              end if
           end do loop_30
           if (d(n) == czero) then
              info = n
              return
           end if
           ! back solve with the matrix u from the factorization.
           do j = 1,nrhs
              b(n,j) = b(n,j)/d(n)
              if (n > 1) b(n - 1,j) = (b(n - 1,j) - du(n - 1)*b(n,j))/d(n - 1)
              do k = n - 2,1,-1
                 b(k,j) = (b(k,j) - du(k)*b(k + 1,j) - dl(k)*b(k + 2,j))/d(k)

              end do
           end do
           return
     end subroutine la_wgtsv
#endif

     !> CGBSV: computes the solution to a complex system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_cgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_cgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_cgbsv
     !> ZGBSV: computes the solution to a complex system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_zgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_zgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_zgbsv
#ifdef LA_WITH_XDP
     !> YGBSV: computes the solution to a complex system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_ygbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_ygbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ygbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_ygbsv
#endif
#ifdef LA_WITH_QP
     !> WGBSV: computes the solution to a complex system of linear equations
     !> A * X = B, where A is a band matrix of order N with KL subdiagonals
     !> and KU superdiagonals, and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as A = L * U, where L is a product of permutation
     !> and unit lower triangular matrices with KL subdiagonals, and U is
     !> upper triangular with KL+KU superdiagonals.  The factored form of A
     !> is then used to solve the system of equations A * X = B.

     pure subroutine la_wgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (kl < 0) then
              info = -2
           else if (ku < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < 2*kl + ku + 1) then
              info = -6
           else if (ldb < max(n,1)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WGBSV ',-info)
              return
           end if
           ! compute the lu factorization of the band matrix a.
           call la_wgbtrf(n,n,kl,ku,ab,ldab,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wgbtrs('NO TRANSPOSE',n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

           end if
           return
     end subroutine la_wgbsv
#endif

     !> CGBSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(inout) :: c(*),r(*)
           complex(sp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================
        ! moved setting of info = n+1 so info does not subsequently get
        ! overwritten.  sven, 17 mar 05.
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(sp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_cgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_claqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_ccopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_cgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_clantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_clangb(norm,n,kl,ku,ab,ldab,rwork)
           rpvgrw = la_clantb('M','U','N',n,kl + ku,afb,ldafb,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_clangb('M',n,kl,ku,ab,ldab,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_cgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,rwork,info)

           ! compute the solution matrix x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_cgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_cgbsvx
     !> ZGBSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(inout) :: c(*),r(*)
           complex(dp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================
        ! moved setting of info = n+1 so info does not subsequently get
        ! overwritten.  sven, 17 mar 05.
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(dp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_zgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_zlaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_zcopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_zgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_zlantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_zlangb(norm,n,kl,ku,ab,ldab,rwork)
           rpvgrw = la_zlantb('M','U','N',n,kl + ku,afb,ldafb,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_zlangb('M',n,kl,ku,ab,ldab,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_zgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,rwork,info)

           ! compute the solution matrix x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_zgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_zgbsvx
#ifdef LA_WITH_XDP
     !> YGBSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_ygbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(xdp),intent(inout) :: c(*),r(*)
           complex(xdp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================
        ! moved setting of info = n+1 so info does not subsequently get
        ! overwritten.  sven, 17 mar 05.
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(xdp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_xlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_ygbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_ylaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_ycopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_ygbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_ylantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_ylangb(norm,n,kl,ku,ab,ldab,rwork)
           rpvgrw = la_ylantb('M','U','N',n,kl + ku,afb,ldafb,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_ylangb('M',n,kl,ku,ab,ldab,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_ygbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,rwork,info)

           ! compute the solution matrix x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ygbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_ygbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_ygbsvx
#endif
#ifdef LA_WITH_QP
     !> WGBSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a band matrix of order N with KL subdiagonals and KU
     !> superdiagonals, and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wgbsvx(fact,trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,equed,r, &
               c,b,ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldafb,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(inout) :: c(*),r(*)
           complex(qp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================
        ! moved setting of info = n+1 so info does not subsequently get
        ! overwritten.  sven, 17 mar 05.
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j,j1,j2
           real(qp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kl < 0) then
              info = -4
           else if (ku < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kl + ku + 1) then
              info = -8
           else if (ldafb < 2*kl + ku + 1) then
              info = -10
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -12
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -13
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -14
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -16
                 else if (ldx < max(1,n)) then
                    info = -18
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_wgbequ(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,infequ)

              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_wlaqgb(n,n,kl,ku,ab,ldab,r,c,rowcnd,colcnd,amax,equed)

                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of the band matrix a.
              do j = 1,n
                 j1 = max(j - ku,1)
                 j2 = min(j + kl,n)
                 call la_wcopy(j2 - j1 + 1,ab(ku + 1 - j + j1,j),1,afb(kl + ku + 1 - j + j1,j),1)

              end do
              call la_wgbtrf(n,n,kl,ku,afb,ldafb,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 anorm = zero
                 do j = 1,info
                    do i = max(ku + 2 - j,1),min(n + ku + 1 - j,kl + ku + 1)
                       anorm = max(anorm,abs(ab(i,j)))
                    end do
                 end do
                 rpvgrw = la_wlantb('M','U','N',info,min(info - 1,kl + ku),afb(max(1, &
                           kl + ku + 2 - info),1),ldafb,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = anorm/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_wlangb(norm,n,kl,ku,ab,ldab,rwork)
           rpvgrw = la_wlantb('M','U','N',n,kl + ku,afb,ldafb,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_wlangb('M',n,kl,ku,ab,ldab,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_wgbcon(norm,n,kl,ku,afb,ldafb,ipiv,anorm,rcond,work,rwork,info)

           ! compute the solution matrix x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wgbtrs(trans,n,kl,ku,nrhs,afb,ldafb,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_wgbrfs(trans,n,kl,ku,nrhs,ab,ldab,afb,ldafb,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_wgbsvx
#endif

     !> CGTSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_cgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           complex(sp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('CGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_ccopy(n,d,1,df,1)
              if (n > 1) then
                 call la_ccopy(n - 1,dl,1,dlf,1)
                 call la_ccopy(n - 1,du,1,duf,1)
              end if
              call la_cgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_clangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_cgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_cgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cgtsvx
     !> ZGTSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_zgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           complex(dp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('ZGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_zcopy(n,d,1,df,1)
              if (n > 1) then
                 call la_zcopy(n - 1,dl,1,dlf,1)
                 call la_zcopy(n - 1,du,1,duf,1)
              end if
              call la_zgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_zlangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_zgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zgtsvx
#ifdef LA_WITH_XDP
     !> YGTSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_ygtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           complex(xdp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('YGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_ycopy(n,d,1,df,1)
              if (n > 1) then
                 call la_ycopy(n - 1,dl,1,dlf,1)
                 call la_ycopy(n - 1,du,1,duf,1)
              end if
              call la_ygttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_ylangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_ygtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ygttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_ygtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_ygtsvx
#endif
#ifdef LA_WITH_QP
     !> WGTSVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     !> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_wgtsvx(fact,trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b, &
               ldb,x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: b(ldb,*),d(*),dl(*),du(*)
           complex(qp),intent(inout) :: df(*),dlf(*),du2(*),duf(*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact,notran
           character :: norm
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           notran = la_lsame(trans,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -14
           else if (ldx < max(1,n)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('WGTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the lu factorization of a.
              call la_wcopy(n,d,1,df,1)
              if (n > 1) then
                 call la_wcopy(n - 1,dl,1,dlf,1)
                 call la_wcopy(n - 1,du,1,duf,1)
              end if
              call la_wgttrf(n,dlf,df,duf,du2,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_wlangt(norm,n,dl,d,du)
           ! compute the reciprocal of the condition number of a.
           call la_wgtcon(norm,n,dlf,df,duf,du2,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wgttrs(trans,n,nrhs,dlf,df,duf,du2,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_wgtrfs(trans,n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,b,ldb,x,ldx, &
                     ferr,berr,work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wgtsvx
#endif

     !> ZCGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> ZCGESV first attempts to factorize the matrix in COMPLEX and use this
     !> factorization within an iterative refinement procedure to produce a
     !> solution with COMPLEX*16 normwise backward error quality (see below).
     !> If the approach fails the method switches to a COMPLEX*16
     !> factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio COMPLEX performance over COMPLEX*16 performance is too
     !> small. A reasonable strategy should take the number of right-hand
     !> sides and the size of the matrix into account. This might be done
     !> with a call to ILAENV in the future. Up to now, we always try
     !> iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by DLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_zcgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,rwork,iter, &
               info)
        use la_constants_dp,only:cone,cnegone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(out) :: rwork(*)
           complex(sp),intent(out) :: swork(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: b(ldb,*)
           complex(dp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(dp),parameter :: bwdmax = 1.0e+00_dp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(dp) :: anrm,cte,eps,rnrm,xnrm
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZCGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip single precision iterative refinement if a priori slower
           ! than double precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_zlange('I',n,n,a,lda,rwork)
           eps = la_dlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=dp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from double precision to single precision and store the
           ! result in sx.
           call la_zlag2c(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from double precision to single precision and store the
           ! result in sa.
           call la_zlag2c(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_cgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_cgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to double precision
           call la_clag2z(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_zlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x,ldx, &
                     cone,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = cabs1(x(la_izamax(n,x(1,i),1),i))
              rnrm = cabs1(work(la_izamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from double precision to single precision
              ! and store the result in sx.
              call la_zlag2c(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_cgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to double precision and update the current
              ! iterate.
              call la_clag2z(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_zaxpy(n,cone,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_zlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x, &
                        ldx,cone,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = cabs1(x(la_izamax(n,x(1,i),1),i))
                 rnrm = cabs1(work(la_izamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the stopping
           ! criterion, set up the iter flag accordingly and follow up on double
           ! precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to double precision.
           call la_zgetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_zlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_zgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_zcgesv
#ifdef LA_WITH_XDP
     !> YZGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> YZGESV first attempts to factorize the matrix in COMPLEX and use this
     !> factorization within an iterative refinement procedure to produce a
     !> solution with COMPLEX*16 normwise backward error quality (see below).
     !> If the approach fails the method switches to a COMPLEX*16
     !> factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio COMPLEX performance over COMPLEX*16 performance is too
     !> small. A reasonable strategy should take the number of right-hand
     !> sides and the size of the matrix into account. This might be done
     !> with a call to ILAENV in the future. Up to now, we always try
     !> iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by XLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_yzgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,rwork,iter, &
               info)
        use la_constants_xdp,only:cone,cnegone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: swork(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: b(ldb,*)
           complex(xdp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(xdp),parameter :: bwdmax = 1.0e+00_xdp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(xdp) :: anrm,cte,eps,rnrm,xnrm
           complex(xdp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YZGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than extended precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_ylange('I',n,n,a,lda,rwork)
           eps = la_xlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=xdp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from extended precision to double precision and store the
           ! result in sx.
           call la_ylag2z(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from extended precision to double precision and store the
           ! result in sa.
           call la_ylag2z(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_zgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_zgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to extended precision
           call la_zlag2y(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_ylacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x,ldx, &
                     cone,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = cabs1(x(la_iyamax(n,x(1,i),1),i))
              rnrm = cabs1(work(la_iyamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from extended precision to double precision
              ! and store the result in sx.
              call la_ylag2z(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_zgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to extended precision and update the current
              ! iterate.
              call la_zlag2y(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_yaxpy(n,cone,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_ylacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x, &
                        ldx,cone,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = cabs1(x(la_iyamax(n,x(1,i),1),i))
                 rnrm = cabs1(work(la_iyamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the stopping
           ! criterion, set up the iter flag accordingly and follow up on double
           ! precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to extended precision.
           call la_ygetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_ylacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_ygetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_yzgesv
#endif
#ifdef LA_WITH_QP
     !> WZGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> WZGESV first attempts to factorize the matrix in COMPLEX and use this
     !> factorization within an iterative refinement procedure to produce a
     !> solution with COMPLEX*16 normwise backward error quality (see below).
     !> If the approach fails the method switches to a COMPLEX*16
     !> factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio COMPLEX performance over COMPLEX*16 performance is too
     !> small. A reasonable strategy should take the number of right-hand
     !> sides and the size of the matrix into account. This might be done
     !> with a call to ILAENV in the future. Up to now, we always try
     !> iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by QLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_wzgesv(n,nrhs,a,lda,ipiv,b,ldb,x,ldx,work,swork,rwork,iter, &
               info)
        use la_constants_qp,only:cone,cnegone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: swork(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: b(ldb,*)
           complex(qp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(qp),parameter :: bwdmax = 1.0e+00_qp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(qp) :: anrm,cte,eps,rnrm,xnrm
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WZGESV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than quad precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_wlange('I',n,n,a,lda,rwork)
           eps = la_qlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=qp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from quad precision to double precision and store the
           ! result in sx.
           call la_wlag2z(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from quad precision to double precision and store the
           ! result in sa.
           call la_wlag2z(n,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the lu factorization of sa.
           call la_zgetrf(n,n,swork(ptsa),n,ipiv,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_zgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx),n, &
                     info)
           ! convert sx back to quad precision
           call la_zlag2w(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_wlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x,ldx, &
                     cone,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = cabs1(x(la_iwamax(n,x(1,i),1),i))
              rnrm = cabs1(work(la_iwamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from quad precision to double precision
              ! and store the result in sx.
              call la_wlag2z(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_zgetrs('NO TRANSPOSE',n,nrhs,swork(ptsa),n,ipiv,swork(ptsx), &
                        n,info)
              ! convert sx back to quad precision and update the current
              ! iterate.
              call la_zlag2w(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_waxpy(n,cone,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_wlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',n,nrhs,n,cnegone,a,lda,x, &
                        ldx,cone,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = cabs1(x(la_iwamax(n,x(1,i),1),i))
                 rnrm = cabs1(work(la_iwamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the stopping
           ! criterion, set up the iter flag accordingly and follow up on double
           ! precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to quad precision.
           call la_wgetrf(n,n,a,lda,ipiv,info)
           if (info /= 0) return
           call la_wlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_wgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,x,ldx,info)
           return
     end subroutine la_wzgesv
#endif

     !> CGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_cgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_cgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_cgesv
     !> ZGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_zgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_zgesv
#ifdef LA_WITH_XDP
     !> YGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_ygesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_ygetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ygetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_ygesv
#endif
#ifdef LA_WITH_QP
     !> WGESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> The LU decomposition with partial pivoting and row interchanges is
     !> used to factor A as
     !> A = P * L * U,
     !> where P is a permutation matrix, L is unit lower triangular, and U is
     !> upper triangular.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_wgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WGESV ',-info)
              return
           end if
           ! compute the lu factorization of a.
           call la_wgetrf(n,n,a,lda,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wgetrs('NO TRANSPOSE',n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_wgesv
#endif

     !> CGESVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(inout) :: c(*),r(*)
           complex(sp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_cgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_claqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_clacpy('FULL',n,n,a,lda,af,ldaf)
              call la_cgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_clantr('M','U','N',info,info,af,ldaf,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_clange('M',n,info,a,lda,rwork)/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_clange(norm,n,n,a,lda,rwork)
           rpvgrw = la_clantr('M','U','N',n,n,af,ldaf,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_clange('M',n,n,a,lda,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_cgecon(norm,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_cgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_cgesvx
     !> ZGESVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(inout) :: c(*),r(*)
           complex(dp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_zgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_zlaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_zlacpy('FULL',n,n,a,lda,af,ldaf)
              call la_zgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_zlantr('M','U','N',info,info,af,ldaf,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_zlange('M',n,info,a,lda,rwork)/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_zlange(norm,n,n,a,lda,rwork)
           rpvgrw = la_zlantr('M','U','N',n,n,af,ldaf,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_zlange('M',n,n,a,lda,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_zgecon(norm,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_zgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_zgesvx
#ifdef LA_WITH_XDP
     !> YGESVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_ygesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(xdp),intent(inout) :: c(*),r(*)
           complex(xdp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(xdp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_xlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_ygeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_ylaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_ylacpy('FULL',n,n,a,lda,af,ldaf)
              call la_ygetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_ylantr('M','U','N',info,info,af,ldaf,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_ylange('M',n,info,a,lda,rwork)/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_ylange(norm,n,n,a,lda,rwork)
           rpvgrw = la_ylantr('M','U','N',n,n,af,ldaf,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_ylange('M',n,n,a,lda,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_ygecon(norm,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ygetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_ygerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_ygesvx
#endif
#ifdef LA_WITH_QP
     !> WGESVX: uses the LU factorization to compute the solution to a complex
     !> system of linear equations
     !> A * X = B,
     !> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wgesvx(fact,trans,n,nrhs,a,lda,af,ldaf,ipiv,equed,r,c,b,ldb, &
               x,ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(inout) :: c(*),r(*)
           complex(qp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: colequ,equil,nofact,notran,rowequ
           character :: norm
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,colcnd,rcmax,rcmin,rowcnd,rpvgrw,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           notran = la_lsame(trans,'N')
           if (nofact .or. equil) then
              equed = 'N'
              rowequ = .false.
              colequ = .false.
           else
              rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
              colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rowequ .or. colequ .or. la_lsame( &
                     equed,'N'))) then
              info = -10
           else
              if (rowequ) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,r(j))
                    rcmax = max(rcmax,r(j))
                 end do
                 if (rcmin <= zero) then
                    info = -11
                 else if (n > 0) then
                    rowcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    rowcnd = one
                 end if
              end if
              if (colequ .and. info == 0) then
                 rcmin = bignum
                 rcmax = zero
                 do j = 1,n
                    rcmin = min(rcmin,c(j))
                    rcmax = max(rcmax,c(j))
                 end do
                 if (rcmin <= zero) then
                    info = -12
                 else if (n > 0) then
                    colcnd = max(rcmin,smlnum)/min(rcmax,bignum)
                 else
                    colcnd = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -14
                 else if (ldx < max(1,n)) then
                    info = -16
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WGESVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_wgeequ(n,n,a,lda,r,c,rowcnd,colcnd,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_wlaqge(n,n,a,lda,r,c,rowcnd,colcnd,amax,equed)
                 rowequ = la_lsame(equed,'R') .or. la_lsame(equed,'B')
                 colequ = la_lsame(equed,'C') .or. la_lsame(equed,'B')
              end if
           end if
           ! scale the right hand side.
           if (notran) then
              if (rowequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       b(i,j) = r(i)*b(i,j)
                    end do
                 end do
              end if
           else if (colequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = c(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the lu factorization of a.
              call la_wlacpy('FULL',n,n,a,lda,af,ldaf)
              call la_wgetrf(n,n,af,ldaf,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 ! compute the reciprocal pivot growth factor of the
                 ! leading rank-deficient info columns of a.
                 rpvgrw = la_wlantr('M','U','N',info,info,af,ldaf,rwork)
                 if (rpvgrw == zero) then
                    rpvgrw = one
                 else
                    rpvgrw = la_wlange('M',n,info,a,lda,rwork)/rpvgrw
                 end if
                 rwork(1) = rpvgrw
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a and the
           ! reciprocal pivot growth factor rpvgrw.
           if (notran) then
              norm = '1'
           else
              norm = 'I'
           end if
           anorm = la_wlange(norm,n,n,a,lda,rwork)
           rpvgrw = la_wlantr('M','U','N',n,n,af,ldaf,rwork)
           if (rpvgrw == zero) then
              rpvgrw = one
           else
              rpvgrw = la_wlange('M',n,n,a,lda,rwork)/rpvgrw
           end if
           ! compute the reciprocal of the condition number of a.
           call la_wgecon(norm,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wgetrs(trans,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_wgerfs(trans,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (notran) then
              if (colequ) then
                 do j = 1,nrhs
                    do i = 1,n
                       x(i,j) = c(i)*x(i,j)
                    end do
                 end do
                 do j = 1,nrhs
                    ferr(j) = ferr(j)/colcnd
                 end do
              end if
           else if (rowequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = r(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/rowcnd
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           rwork(1) = rpvgrw
           return
     end subroutine la_wgesvx
#endif

end module la_lapack_solve_lu
