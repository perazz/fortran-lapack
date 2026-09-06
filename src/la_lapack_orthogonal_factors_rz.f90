!> RZ factorization: trapezoidal reduction and its reflectors
module la_lapack_orthogonal_factors_rz
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_blas_like_l1
     use la_lapack_householder_reflectors
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slarz
     public :: la_slarzb
     public :: la_slarzt
     public :: la_sormr3
     public :: la_sormrz
     public :: la_slatrz
     public :: la_stzrzf
     public :: la_dlarz
     public :: la_dlarzb
     public :: la_dlarzt
     public :: la_dormr3
     public :: la_dormrz
     public :: la_dlatrz
     public :: la_dtzrzf
     public :: la_qlarz
     public :: la_qlarzb
     public :: la_qlarzt
     public :: la_qormr3
     public :: la_qormrz
     public :: la_qlatrz
     public :: la_qtzrzf
     public :: la_clarz
     public :: la_clarzb
     public :: la_clarzt
     public :: la_clatrz
     public :: la_ctzrzf
     public :: la_cunmr3
     public :: la_cunmrz
     public :: la_zlarz
     public :: la_zlarzb
     public :: la_zlarzt
     public :: la_zlatrz
     public :: la_ztzrzf
     public :: la_zunmr3
     public :: la_zunmrz
     public :: la_wlarz
     public :: la_wlarzb
     public :: la_wlarzt
     public :: la_wlatrz
     public :: la_wtzrzf
     public :: la_wunmr3
     public :: la_wunmrz

     contains

     !> SLARZ: applies a real elementary reflector H to a real M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> H is a product of k elementary reflectors as returned by STZRZF.

     pure subroutine la_slarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           real(sp),intent(in) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(in) :: v(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= zero) then
                 ! w( 1:n ) = c( 1, 1:n )
                 call la_scopy(n,c,ldc,work,1)
                 ! w( 1:n ) = w( 1:n ) + c( m-l+1:m, 1:n )**t * v( 1:l )
                 call la_sgemv('TRANSPOSE',l,n,one,c(m - l + 1,1),ldc,v,incv,one,work, &
                            1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_saxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**t
                 call la_sger(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= zero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_scopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_sgemv('NO TRANSPOSE',m,l,one,c(1,n - l + 1),ldc,v,incv,one, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_saxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**t
                 call la_sger(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_slarz
     !> DLARZ: applies a real elementary reflector H to a real M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> H is a product of k elementary reflectors as returned by DTZRZF.

     pure subroutine la_dlarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           real(dp),intent(in) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(in) :: v(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= zero) then
                 ! w( 1:n ) = c( 1, 1:n )
                 call la_dcopy(n,c,ldc,work,1)
                 ! w( 1:n ) = w( 1:n ) + c( m-l+1:m, 1:n )**t * v( 1:l )
                 call la_dgemv('TRANSPOSE',l,n,one,c(m - l + 1,1),ldc,v,incv,one,work, &
                            1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_daxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**t
                 call la_dger(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= zero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_dcopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_dgemv('NO TRANSPOSE',m,l,one,c(1,n - l + 1),ldc,v,incv,one, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_daxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**t
                 call la_dger(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_dlarz
     !> QLARZ: applies a real elementary reflector H to a real M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> H is a product of k elementary reflectors as returned by QTZRZF.

     pure subroutine la_qlarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           real(qp),intent(in) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(in) :: v(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= zero) then
                 ! w( 1:n ) = c( 1, 1:n )
                 call la_qcopy(n,c,ldc,work,1)
                 ! w( 1:n ) = w( 1:n ) + c( m-l+1:m, 1:n )**t * v( 1:l )
                 call la_qgemv('TRANSPOSE',l,n,one,c(m - l + 1,1),ldc,v,incv,one,work, &
                            1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_qaxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**t
                 call la_qger(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= zero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_qcopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_qgemv('NO TRANSPOSE',m,l,one,c(1,n - l + 1),ldc,v,incv,one, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_qaxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**t
                 call la_qger(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_qlarz

     !> SLARZB: applies a real block reflector H or its transpose H**T to
     !> a real distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_slarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           real(sp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**t * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**t
              do j = 1,k
                 call la_scopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**t * v( 1:k, 1:l )**t
              if (l > 0) call la_sgemm('TRANSPOSE','TRANSPOSE',n,k,l,one,c(m - l + 1,1), &
                        ldc,v,ldv,one,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_strmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**t
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**t * w( 1:n, 1:k )**t
              if (l > 0) call la_sgemm('TRANSPOSE','TRANSPOSE',l,n,k,-one,v,ldv,work, &
                        ldwork,one,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**t
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_scopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**t
              if (l > 0) call la_sgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,one,c(1,n - l + 1), &
                         ldc,v,ldv,one,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * t  or  w( 1:m, 1:k ) * t**t
              call la_strmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * v( 1:k, 1:l )
              if (l > 0) call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-one,work, &
                        ldwork,v,ldv,one,c(1,n - l + 1),ldc)
           end if
           return
     end subroutine la_slarzb
     !> DLARZB: applies a real block reflector H or its transpose H**T to
     !> a real distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_dlarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           real(dp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**t * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**t
              do j = 1,k
                 call la_dcopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**t * v( 1:k, 1:l )**t
              if (l > 0) call la_dgemm('TRANSPOSE','TRANSPOSE',n,k,l,one,c(m - l + 1,1), &
                        ldc,v,ldv,one,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_dtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**t
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**t * w( 1:n, 1:k )**t
              if (l > 0) call la_dgemm('TRANSPOSE','TRANSPOSE',l,n,k,-one,v,ldv,work, &
                        ldwork,one,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**t
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_dcopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**t
              if (l > 0) call la_dgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,one,c(1,n - l + 1), &
                         ldc,v,ldv,one,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * t  or  w( 1:m, 1:k ) * t**t
              call la_dtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * v( 1:k, 1:l )
              if (l > 0) call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-one,work, &
                        ldwork,v,ldv,one,c(1,n - l + 1),ldc)
           end if
           return
     end subroutine la_dlarzb
     !> QLARZB: applies a real block reflector H or its transpose H**T to
     !> a real distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_qlarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           real(qp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**t * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**t
              do j = 1,k
                 call la_qcopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**t * v( 1:k, 1:l )**t
              if (l > 0) call la_qgemm('TRANSPOSE','TRANSPOSE',n,k,l,one,c(m - l + 1,1), &
                        ldc,v,ldv,one,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_qtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**t
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**t * w( 1:n, 1:k )**t
              if (l > 0) call la_qgemm('TRANSPOSE','TRANSPOSE',l,n,k,-one,v,ldv,work, &
                        ldwork,one,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**t
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_qcopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**t
              if (l > 0) call la_qgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,one,c(1,n - l + 1), &
                         ldc,v,ldv,one,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * t  or  w( 1:m, 1:k ) * t**t
              call la_qtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt,work, &
                        ldwork)
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * v( 1:k, 1:l )
              if (l > 0) call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-one,work, &
                        ldwork,v,ldv,one,c(1,n - l + 1),ldc)
           end if
           return
     end subroutine la_qlarzb

     !> SLARZT: forms the triangular factor T of a real block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_slarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(sp),intent(out) :: t(ldt,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('SLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == zero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = zero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**t
                    call la_sgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,zero,t(i + 1,i),1)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_strmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_slarzt
     !> DLARZT: forms the triangular factor T of a real block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_dlarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(dp),intent(out) :: t(ldt,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('DLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == zero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = zero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**t
                    call la_dgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,zero,t(i + 1,i),1)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_dtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_dlarzt
     !> QLARZT: forms the triangular factor T of a real block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_qlarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(qp),intent(out) :: t(ldt,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('QLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == zero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = zero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**t
                    call la_qgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,zero,t(i + 1,i),1)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_qtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_qlarzt

     !> SORMR3: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'C',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by STZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_sormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),tau(*)
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('SORMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**t is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**t is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**t
              call la_slarz(side,mi,ni,l,a(i,ja),lda,tau(i),c(ic,jc),ldc, &
                        work)
           end do
           return
     end subroutine la_sormr3
     !> DORMR3: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'C',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by DTZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_dormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),tau(*)
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DORMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**t is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**t is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**t
              call la_dlarz(side,mi,ni,l,a(i,ja),lda,tau(i),c(ic,jc),ldc, &
                        work)
           end do
           return
     end subroutine la_dormr3
     !> QORMR3: overwrites the general real m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**T* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**T if SIDE = 'R' and TRANS = 'C',
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by QTZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_qormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),tau(*)
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QORMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
              ic = 1
           end if
           do i = i1,i2,i3
              if (left) then
                 ! h(i) or h(i)**t is applied to c(i:m,1:n)
                 mi = m - i + 1
                 ic = i
              else
                 ! h(i) or h(i)**t is applied to c(1:m,i:n)
                 ni = n - i + 1
                 jc = i
              end if
              ! apply h(i) or h(i)**t
              call la_qlarz(side,mi,ni,l,a(i,ja),lda,tau(i),c(ic,jc),ldc, &
                        work)
           end do
           return
     end subroutine la_qormr3

     !> SORMRZ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by STZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_sormrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'SORMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORMRZ',-info)
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
                 nbmin = max(2,la_ilaenv(2,'SORMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_sormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_slarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_slarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_sormrz
     !> DORMRZ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_dormrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'DORMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORMRZ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'DORMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_dormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_dlarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_dlarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_dormrz
     !> QORMRZ: overwrites the general real M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> where Q is a real orthogonal matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by QTZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_qormrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'QORMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORMRZ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'QORMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_qormr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_qlarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_qlarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_qormrz

     !> SLATRZ: factors the M-by-(M+L) real upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means
     !> of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_slatrz(m,n,l,a,lda,tau,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! test the input arguments
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_slarfg(l + 1,a(i,i),a(i,n - l + 1),lda,tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_slarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,tau(i),a(1,i), &
                        lda,work)
           end do
           return
     end subroutine la_slatrz
     !> DLATRZ: factors the M-by-(M+L) real upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means
     !> of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_dlatrz(m,n,l,a,lda,tau,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! test the input arguments
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_dlarfg(l + 1,a(i,i),a(i,n - l + 1),lda,tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_dlarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,tau(i),a(1,i), &
                        lda,work)
           end do
           return
     end subroutine la_dlatrz
     !> QLATRZ: factors the M-by-(M+L) real upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z, by means
     !> of orthogonal transformations.  Z is an (M+L)-by-(M+L) orthogonal
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_qlatrz(m,n,l,a,lda,tau,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Executable Statements
           ! test the input arguments
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_qlarfg(l + 1,a(i,i),a(i,n - l + 1),lda,tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_qlarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,tau(i),a(1,i), &
                        lda,work)
           end do
           return
     end subroutine la_qlatrz

     !> STZRZF: reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
     !> to upper triangular form by means of orthogonal transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_stzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_sp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'SGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('STZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'SGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'SGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_slatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_slarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_slarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_slatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_stzrzf
     !> DTZRZF: reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
     !> to upper triangular form by means of orthogonal transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_dtzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_dp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'DGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DTZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'DGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'DGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_dlatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_dlarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_dlarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_dlatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_dtzrzf
     !> QTZRZF: reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
     !> to upper triangular form by means of orthogonal transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_qtzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_qp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'QGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QTZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = zero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'QGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'QGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_qlatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_qlarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_qlarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_qlatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_qtzrzf

     !> CLARZ: applies a complex elementary reflector H to a complex
     !> M-by-N matrix C, from either the left or the right. H is represented
     !> in the form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
     !> tau.
     !> H is a product of k elementary reflectors as returned by CTZRZF.

     pure subroutine la_clarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           complex(sp),intent(in) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(in) :: v(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= czero) then
                 ! w( 1:n ) = conjg( c( 1, 1:n ) )
                 call la_ccopy(n,c,ldc,work,1)
                 call la_clacgv(n,work,1)
                 ! w( 1:n ) = conjg( w( 1:n ) + c( m-l+1:m, 1:n )**h * v( 1:l ) )
                 call la_cgemv('CONJUGATE TRANSPOSE',l,n,cone,c(m - l + 1,1),ldc,v,incv, &
                            cone,work,1)
                 call la_clacgv(n,work,1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_caxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**h
                 call la_cgeru(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= czero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_ccopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_cgemv('NO TRANSPOSE',m,l,cone,c(1,n - l + 1),ldc,v,incv,cone, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_caxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**h
                 call la_cgerc(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_clarz
     !> ZLARZ: applies a complex elementary reflector H to a complex
     !> M-by-N matrix C, from either the left or the right. H is represented
     !> in the form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
     !> tau.
     !> H is a product of k elementary reflectors as returned by ZTZRZF.

     pure subroutine la_zlarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           complex(dp),intent(in) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(in) :: v(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= czero) then
                 ! w( 1:n ) = conjg( c( 1, 1:n ) )
                 call la_zcopy(n,c,ldc,work,1)
                 call la_zlacgv(n,work,1)
                 ! w( 1:n ) = conjg( w( 1:n ) + c( m-l+1:m, 1:n )**h * v( 1:l ) )
                 call la_zgemv('CONJUGATE TRANSPOSE',l,n,cone,c(m - l + 1,1),ldc,v,incv, &
                            cone,work,1)
                 call la_zlacgv(n,work,1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_zaxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**h
                 call la_zgeru(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= czero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_zcopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_zgemv('NO TRANSPOSE',m,l,cone,c(1,n - l + 1),ldc,v,incv,cone, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_zaxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**h
                 call la_zgerc(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_zlarz
     !> WLARZ: applies a complex elementary reflector H to a complex
     !> M-by-N matrix C, from either the left or the right. H is represented
     !> in the form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
     !> tau.
     !> H is a product of k elementary reflectors as returned by WTZRZF.

     pure subroutine la_wlarz(side,m,n,l,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,l,ldc,m,n
           complex(qp),intent(in) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(in) :: v(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Executable Statements
           if (la_lsame(side,'L')) then
              ! form  h * c
              if (tau /= czero) then
                 ! w( 1:n ) = conjg( c( 1, 1:n ) )
                 call la_wcopy(n,c,ldc,work,1)
                 call la_wlacgv(n,work,1)
                 ! w( 1:n ) = conjg( w( 1:n ) + c( m-l+1:m, 1:n )**h * v( 1:l ) )
                 call la_wgemv('CONJUGATE TRANSPOSE',l,n,cone,c(m - l + 1,1),ldc,v,incv, &
                            cone,work,1)
                 call la_wlacgv(n,work,1)
                 ! c( 1, 1:n ) = c( 1, 1:n ) - tau * w( 1:n )
                 call la_waxpy(n,-tau,work,1,c,ldc)
                 ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                     ! tau * v( 1:l ) * w( 1:n )**h
                 call la_wgeru(l,n,-tau,v,incv,work,1,c(m - l + 1,1),ldc)
              end if
           else
              ! form  c * h
              if (tau /= czero) then
                 ! w( 1:m ) = c( 1:m, 1 )
                 call la_wcopy(m,c,1,work,1)
                 ! w( 1:m ) = w( 1:m ) + c( 1:m, n-l+1:n, 1:n ) * v( 1:l )
                 call la_wgemv('NO TRANSPOSE',m,l,cone,c(1,n - l + 1),ldc,v,incv,cone, &
                           work,1)
                 ! c( 1:m, 1 ) = c( 1:m, 1 ) - tau * w( 1:m )
                 call la_waxpy(m,-tau,work,1,c,1)
                 ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                     ! tau * w( 1:m ) * v( 1:l )**h
                 call la_wgerc(m,l,-tau,work,1,v,incv,c(1,n - l + 1),ldc)
              end if
           end if
           return
     end subroutine la_wlarz

     !> CLARZB: applies a complex block reflector H or its transpose H**H
     !> to a complex distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_clarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           complex(sp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**h * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**h
              do j = 1,k
                 call la_ccopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**h * v( 1:k, 1:l )**t
              if (l > 0) call la_cgemm('TRANSPOSE','CONJUGATE TRANSPOSE',n,k,l,cone,c(m - &
                        l + 1,1),ldc,v,ldv,cone,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_ctrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**h
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**h * w( 1:n, 1:k )**h
              if (l > 0) call la_cgemm('TRANSPOSE','TRANSPOSE',l,n,k,-cone,v,ldv,work, &
                        ldwork,cone,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**h
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_ccopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**h
              if (l > 0) call la_cgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,cone,c(1,n - l + 1) &
                        ,ldc,v,ldv,cone,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * conjg( t )  or
                              ! w( 1:m, 1:k ) * t**h
              do j = 1,k
                 call la_clacgv(k - j + 1,t(j,j),1)
              end do
              call la_ctrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt,work, &
                        ldwork)
              do j = 1,k
                 call la_clacgv(k - j + 1,t(j,j),1)
              end do
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * conjg( v( 1:k, 1:l ) )
              do j = 1,l
                 call la_clacgv(k,v(1,j),1)
              end do
              if (l > 0) call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-cone,work, &
                        ldwork,v,ldv,cone,c(1,n - l + 1),ldc)
              do j = 1,l
                 call la_clacgv(k,v(1,j),1)
              end do
           end if
           return
     end subroutine la_clarzb
     !> ZLARZB: applies a complex block reflector H or its transpose H**H
     !> to a complex distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_zlarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           complex(dp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**h * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**h
              do j = 1,k
                 call la_zcopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**h * v( 1:k, 1:l )**t
              if (l > 0) call la_zgemm('TRANSPOSE','CONJUGATE TRANSPOSE',n,k,l,cone,c(m - &
                        l + 1,1),ldc,v,ldv,cone,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_ztrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**h
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**h * w( 1:n, 1:k )**h
              if (l > 0) call la_zgemm('TRANSPOSE','TRANSPOSE',l,n,k,-cone,v,ldv,work, &
                        ldwork,cone,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**h
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_zcopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**h
              if (l > 0) call la_zgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,cone,c(1,n - l + 1) &
                        ,ldc,v,ldv,cone,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * conjg( t )  or
                              ! w( 1:m, 1:k ) * t**h
              do j = 1,k
                 call la_zlacgv(k - j + 1,t(j,j),1)
              end do
              call la_ztrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt,work, &
                        ldwork)
              do j = 1,k
                 call la_zlacgv(k - j + 1,t(j,j),1)
              end do
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * conjg( v( 1:k, 1:l ) )
              do j = 1,l
                 call la_zlacgv(k,v(1,j),1)
              end do
              if (l > 0) call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-cone,work, &
                        ldwork,v,ldv,cone,c(1,n - l + 1),ldc)
              do j = 1,l
                 call la_zlacgv(k,v(1,j),1)
              end do
           end if
           return
     end subroutine la_zlarzb
     !> WLARZB: applies a complex block reflector H or its transpose H**H
     !> to a complex distributed M-by-N  C from the left or the right.
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_wlarzb(side,trans,direct,storev,m,n,k,l,v,ldv,t,ldt,c, &
               ldc,work,ldwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,l,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*),t(ldt,*),v(ldv,*)
           complex(qp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -3
           else if (.not. la_lsame(storev,'R')) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WLARZB',-info)
              return
           end if
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(side,'L')) then
              ! form  h * c  or  h**h * c
              ! w( 1:n, 1:k ) = c( 1:k, 1:n )**h
              do j = 1,k
                 call la_wcopy(n,c(j,1),ldc,work(1,j),1)
              end do
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) + ...
                              ! c( m-l+1:m, 1:n )**h * v( 1:k, 1:l )**t
              if (l > 0) call la_wgemm('TRANSPOSE','CONJUGATE TRANSPOSE',n,k,l,cone,c(m - &
                        l + 1,1),ldc,v,ldv,cone,work,ldwork)
              ! w( 1:n, 1:k ) = w( 1:n, 1:k ) * t**t  or  w( 1:m, 1:k ) * t
              call la_wtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt,work, &
                        ldwork)
              ! c( 1:k, 1:n ) = c( 1:k, 1:n ) - w( 1:n, 1:k )**h
              do j = 1,n
                 do i = 1,k
                    c(i,j) = c(i,j) - work(j,i)
                 end do
              end do
              ! c( m-l+1:m, 1:n ) = c( m-l+1:m, 1:n ) - ...
                                  ! v( 1:k, 1:l )**h * w( 1:n, 1:k )**h
              if (l > 0) call la_wgemm('TRANSPOSE','TRANSPOSE',l,n,k,-cone,v,ldv,work, &
                        ldwork,cone,c(m - l + 1,1),ldc)
           else if (la_lsame(side,'R')) then
              ! form  c * h  or  c * h**h
              ! w( 1:m, 1:k ) = c( 1:m, 1:k )
              do j = 1,k
                 call la_wcopy(m,c(1,j),1,work(1,j),1)
              end do
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) + ...
                              ! c( 1:m, n-l+1:n ) * v( 1:k, 1:l )**h
              if (l > 0) call la_wgemm('NO TRANSPOSE','TRANSPOSE',m,k,l,cone,c(1,n - l + 1) &
                        ,ldc,v,ldv,cone,work,ldwork)
              ! w( 1:m, 1:k ) = w( 1:m, 1:k ) * conjg( t )  or
                              ! w( 1:m, 1:k ) * t**h
              do j = 1,k
                 call la_wlacgv(k - j + 1,t(j,j),1)
              end do
              call la_wtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt,work, &
                        ldwork)
              do j = 1,k
                 call la_wlacgv(k - j + 1,t(j,j),1)
              end do
              ! c( 1:m, 1:k ) = c( 1:m, 1:k ) - w( 1:m, 1:k )
              do j = 1,k
                 do i = 1,m
                    c(i,j) = c(i,j) - work(i,j)
                 end do
              end do
              ! c( 1:m, n-l+1:n ) = c( 1:m, n-l+1:n ) - ...
                                  ! w( 1:m, 1:k ) * conjg( v( 1:k, 1:l ) )
              do j = 1,l
                 call la_wlacgv(k,v(1,j),1)
              end do
              if (l > 0) call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m,l,k,-cone,work, &
                        ldwork,v,ldv,cone,c(1,n - l + 1),ldc)
              do j = 1,l
                 call la_wlacgv(k,v(1,j),1)
              end do
           end if
           return
     end subroutine la_wlarzb

     !> CLARZT: forms the triangular factor T of a complex block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_clarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(sp),intent(out) :: t(ldt,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('CLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == czero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = czero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**h
                    call la_clacgv(n,v(i,1),ldv)
                    call la_cgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,czero,t(i + 1,i),1)
                    call la_clacgv(n,v(i,1),ldv)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_ctrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_clarzt
     !> ZLARZT: forms the triangular factor T of a complex block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_zlarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(dp),intent(out) :: t(ldt,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('ZLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == czero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = czero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**h
                    call la_zlacgv(n,v(i,1),ldv)
                    call la_zgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,czero,t(i + 1,i),1)
                    call la_zlacgv(n,v(i,1),ldv)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_ztrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_zlarzt
     !> WLARZT: forms the triangular factor T of a complex block reflector
     !> H of order > n, which is defined as a product of k elementary
     !> reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V
     !> Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

     pure subroutine la_wlarzt(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(qp),intent(out) :: t(ldt,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(inout) :: v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           ! Executable Statements
           ! check for currently supported options
           info = 0
           if (.not. la_lsame(direct,'B')) then
              info = -1
           else if (.not. la_lsame(storev,'R')) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('WLARZT',-info)
              return
           end if
           do i = k,1,-1
              if (tau(i) == czero) then
                 ! h(i)  =  i
                 do j = i,k
                    t(j,i) = czero
                 end do
              else
                 ! general case
                 if (i < k) then
                    ! t(i+1:k,i) = - tau(i) * v(i+1:k,1:n) * v(i,1:n)**h
                    call la_wlacgv(n,v(i,1),ldv)
                    call la_wgemv('NO TRANSPOSE',k - i,n,-tau(i),v(i + 1,1),ldv,v(i, &
                              1),ldv,czero,t(i + 1,i),1)
                    call la_wlacgv(n,v(i,1),ldv)
                    ! t(i+1:k,i) = t(i+1:k,i+1:k) * t(i+1:k,i)
                    call la_wtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                              ldt,t(i + 1,i),1)
                 end if
                 t(i,i) = tau(i)
              end if
           end do
           return
     end subroutine la_wlarzt

     !> CLATRZ: factors the M-by-(M+L) complex upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means
     !> of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_clatrz(m,n,l,a,lda,tau,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_clacgv(l,a(i,n - l + 1),lda)
              alpha = conjg(a(i,i))
              call la_clarfg(l + 1,alpha,a(i,n - l + 1),lda,tau(i))
              tau(i) = conjg(tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_clarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,conjg(tau(i)),a( &
                        1,i),lda,work)
              a(i,i) = conjg(alpha)
           end do
           return
     end subroutine la_clatrz
     !> ZLATRZ: factors the M-by-(M+L) complex upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means
     !> of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_zlatrz(m,n,l,a,lda,tau,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_zlacgv(l,a(i,n - l + 1),lda)
              alpha = conjg(a(i,i))
              call la_zlarfg(l + 1,alpha,a(i,n - l + 1),lda,tau(i))
              tau(i) = conjg(tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_zlarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,conjg(tau(i)),a( &
                        1,i),lda,work)
              a(i,i) = conjg(alpha)
           end do
           return
     end subroutine la_zlatrz
     !> WLATRZ: factors the M-by-(M+L) complex upper trapezoidal matrix
     !> [ A1 A2 ] = [ A(1:M,1:M) A(1:M,N-L+1:N) ] as ( R  0 ) * Z by means
     !> of unitary transformations, where  Z is an (M+L)-by-(M+L) unitary
     !> matrix and, R and A1 are M-by-M upper triangular matrices.

     pure subroutine la_wlatrz(m,n,l,a,lda,tau,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: l,lda,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: tau(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           do i = m,1,-1
              ! generate elementary reflector h(i) to annihilate
              ! [ a(i,i) a(i,n-l+1:n) ]
              call la_wlacgv(l,a(i,n - l + 1),lda)
              alpha = conjg(a(i,i))
              call la_wlarfg(l + 1,alpha,a(i,n - l + 1),lda,tau(i))
              tau(i) = conjg(tau(i))
              ! apply h(i) to a(1:i-1,i:n) from the right
              call la_wlarz('RIGHT',i - 1,n - i + 1,l,a(i,n - l + 1),lda,conjg(tau(i)),a( &
                        1,i),lda,work)
              a(i,i) = conjg(alpha)
           end do
           return
     end subroutine la_wlatrz

     !> CTZRZF: reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A
     !> to upper triangular form by means of unitary transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N unitary matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_ctzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_sp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'CGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CTZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'CGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'CGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_clatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_clarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_clarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_clatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_ctzrzf
     !> ZTZRZF: reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A
     !> to upper triangular form by means of unitary transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N unitary matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_ztzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_dp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'ZGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZTZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'ZGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'ZGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_zlatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_zlarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_zlarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_zlatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_ztzrzf
     !> WTZRZF: reduces the M-by-N ( M<=N ) complex upper trapezoidal matrix A
     !> to upper triangular form by means of unitary transformations.
     !> The upper trapezoidal matrix A is factored as
     !> A = ( R  0 ) * Z,
     !> where Z is an N-by-N unitary matrix and R is an M-by-M upper
     !> triangular matrix.

     pure subroutine la_wtzrzf(m,n,a,lda,tau,work,lwork,info)
        use la_constants_qp
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
           integer(ilp) :: i,ib,iws,ki,kk,ldwork,lwkmin,lwkopt,m1,mu,nb,nbmin, &
                     nx
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < m) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info == 0) then
              if (m == 0 .or. m == n) then
                 lwkopt = 1
                 lwkmin = 1
              else
                 ! determine the block size.
                 nb = la_ilaenv(1,'WGERQF',' ',m,n,-1,-1)
                 lwkopt = m*nb
                 lwkmin = max(1,m)
              end if
              work(1) = lwkopt
              if (lwork < lwkmin .and. .not. lquery) then
                 info = -7
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WTZRZF',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0) then
              return
           else if (m == n) then
              do i = 1,n
                 tau(i) = czero
              end do
              return
           end if
           nbmin = 2
           nx = 1
           iws = m
           if (nb > 1 .and. nb < m) then
              ! determine when to cross over from blocked to unblocked code.
              nx = max(0,la_ilaenv(3,'WGERQF',' ',m,n,-1,-1))
              if (nx < m) then
                 ! determine if workspace is large enough for blocked code.
                 ldwork = m
                 iws = ldwork*nb
                 if (lwork < iws) then
                    ! not enough workspace to use optimal nb:  reduce nb and
                    ! determine the minimum value of nb.
                    nb = lwork/ldwork
                    nbmin = max(2,la_ilaenv(2,'WGERQF',' ',m,n,-1,-1))
                 end if
              end if
           end if
           if (nb >= nbmin .and. nb < m .and. nx < m) then
              ! use blocked code initially.
              ! the last kk rows are handled by the block method.
              m1 = min(m + 1,n)
              ki = ((m - nx - 1)/nb)*nb
              kk = min(m,ki + nb)
              do i = m - kk + ki + 1,m - kk + 1,-nb
                 ib = min(m - i + 1,nb)
                 ! compute the tz factorization of the current block
                 ! a(i:i+ib-1,i:n)
                 call la_wlatrz(ib,n - i + 1,n - m,a(i,i),lda,tau(i),work)
                 if (i > 1) then
                    ! form the triangular factor of the block reflector
                    ! h = h(i+ib-1) . . . h(i+1) h(i)
                    call la_wlarzt('BACKWARD','ROWWISE',n - m,ib,a(i,m1),lda,tau(i), &
                              work,ldwork)
                    ! apply h to a(1:i-1,i:n) from the right
                    call la_wlarzb('RIGHT','NO TRANSPOSE','BACKWARD','ROWWISE',i - 1,n - i + 1, &
                     ib,n - m,a(i,m1),lda,work,ldwork,a(1,i),lda,work(ib + 1),ldwork)

                 end if
              end do
              mu = i + nb - 1
           else
              mu = m
           end if
           ! use unblocked code to factor the last or only block
           if (mu > 0) call la_wlatrz(mu,n,n - m,a,lda,tau,work)
           work(1) = lwkopt
           return
     end subroutine la_wtzrzf

     !> CUNMR3: overwrites the general complex m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by CTZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_cunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),tau(*)
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
           complex(sp) :: taui
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CUNMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
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
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              call la_clarz(side,mi,ni,l,a(i,ja),lda,taui,c(ic,jc),ldc,work)

           end do
           return
     end subroutine la_cunmr3
     !> ZUNMR3: overwrites the general complex m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by ZTZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_zunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),tau(*)
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
           complex(dp) :: taui
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
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
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              call la_zlarz(side,mi,ni,l,a(i,ja),lda,taui,c(ic,jc),ldc,work)

           end do
           return
     end subroutine la_zunmr3
     !> WUNMR3: overwrites the general complex m by n matrix C with
     !> Q * C  if SIDE = 'L' and TRANS = 'N', or
     !> Q**H* C  if SIDE = 'L' and TRANS = 'C', or
     !> C * Q  if SIDE = 'R' and TRANS = 'N', or
     !> C * Q**H if SIDE = 'R' and TRANS = 'C',
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by WTZRZF. Q is of order m if SIDE = 'L' and of order n
     !> if SIDE = 'R'.

     pure subroutine la_wunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,m,n
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),tau(*)
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: left,notran
           integer(ilp) :: i,i1,i2,i3,ic,ja,jc,mi,ni,nq
           complex(qp) :: taui
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WUNMR3',-info)
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0 .or. k == 0) return
           if ((left .and. .not. notran .or. .not. left .and. notran)) then
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
              ja = m - l + 1
              jc = 1
           else
              mi = m
              ja = n - l + 1
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
                 taui = tau(i)
              else
                 taui = conjg(tau(i))
              end if
              call la_wlarz(side,mi,ni,l,a(i,ja),lda,taui,c(ic,jc),ldc,work)

           end do
           return
     end subroutine la_wunmr3

     !> CUNMRZ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by CTZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_cunmrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'CUNMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNMRZ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           ! determine the block size.
           nb = min(nbmax,la_ilaenv(1,'CUNMRQ',side//trans,m,n,k,-1))
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'CUNMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_cunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_clarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_clarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_cunmrz
     !> ZUNMRZ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by ZTZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_zunmrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < max(1,nw) .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'ZUNMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMRZ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           ! determine the block size.  nb may be at most nbmax, where nbmax
           ! is used to define the local array t.
           nb = min(nbmax,la_ilaenv(1,'ZUNMRQ',side//trans,m,n,k,-1))
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'ZUNMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_zunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_zlarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_zlarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_zunmrz
     !> WUNMRZ: overwrites the general complex M-by-N matrix C with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> where Q is a complex unitary matrix defined as the product of k
     !> elementary reflectors
     !> Q = H(1) H(2) . . . H(k)
     !> as returned by WTZRZF. Q is of order M if SIDE = 'L' and of order N
     !> if SIDE = 'R'.

     pure subroutine la_wunmrz(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,l,lda,ldc,lwork,m,n
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
           integer(ilp) :: i,i1,i2,i3,ib,ic,iinfo,iwt,ja,jc,ldwork,lwkopt,mi,nb, &
                     nbmin,ni,nq,nw
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
           else if (l < 0 .or. (left .and. (l > m)) .or. (.not. left .and. (l > n))) then
              info = -6
           else if (lda < max(1,k)) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < max(1,nw) .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              ! compute the workspace requirements
              if (m == 0 .or. n == 0) then
                 lwkopt = 1
              else
                 nb = min(nbmax,la_ilaenv(1,'WUNMRQ',side//trans,m,n,k,-1))

                 lwkopt = nw*nb + tsize
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNMRZ',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              return
           end if
           ! determine the block size.  nb may be at most nbmax, where nbmax
           ! is used to define the local array t.
           nb = min(nbmax,la_ilaenv(1,'WUNMRQ',side//trans,m,n,k,-1))
           nbmin = 2
           ldwork = nw
           if (nb > 1 .and. nb < k) then
              if (lwork < lwkopt) then
                 nb = (lwork - tsize)/ldwork
                 nbmin = max(2,la_ilaenv(2,'WUNMRQ',side//trans,m,n,k,-1))
              end if
           end if
           if (nb < nbmin .or. nb >= k) then
              ! use unblocked code
              call la_wunmr3(side,trans,m,n,k,l,a,lda,tau,c,ldc,work,iinfo)

           else
              ! use blocked code
              iwt = 1 + nw*nb
              if ((left .and. .not. notran) .or. (.not. left .and. notran)) then
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
                 ja = m - l + 1
              else
                 mi = m
                 ic = 1
                 ja = n - l + 1
              end if
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              do i = i1,i2,i3
                 ib = min(nb,k - i + 1)
                 ! form the triangular factor of the block reflector
                 ! h = h(i+ib-1) . . . h(i+1) h(i)
                 call la_wlarzt('BACKWARD','ROWWISE',l,ib,a(i,ja),lda,tau(i),work( &
                            iwt),ldt)
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
                 call la_wlarzb(side,transt,'BACKWARD','ROWWISE',mi,ni,ib,l,a(i,ja) &
                           ,lda,work(iwt),ldt,c(ic,jc),ldc,work,ldwork)
              end do
           end if
           work(1) = lwkopt
           return
     end subroutine la_wunmrz

end module la_lapack_orthogonal_factors_rz
