!> Householder reflectors: generation, blocking, application
module la_lapack_householder_reflectors
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_sym
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_scalar
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slarf
     public :: la_slarfb
     public :: la_slarft
     public :: la_slarfx
     public :: la_slarfy
     public :: la_slarfg
     public :: la_slarfgp
     public :: la_dlarf
     public :: la_dlarfb
     public :: la_dlarft
     public :: la_dlarfx
     public :: la_dlarfy
     public :: la_dlarfg
     public :: la_dlarfgp
#ifdef LA_WITH_XDP
     public :: la_xlarf
     public :: la_xlarfb
     public :: la_xlarft
     public :: la_xlarfx
     public :: la_xlarfy
     public :: la_xlarfg
     public :: la_xlarfgp
#endif
#ifdef LA_WITH_QP
     public :: la_qlarf
     public :: la_qlarfb
     public :: la_qlarft
     public :: la_qlarfx
     public :: la_qlarfy
     public :: la_qlarfg
     public :: la_qlarfgp
#endif
     public :: la_clarf
     public :: la_clarfb
     public :: la_clarfg
     public :: la_clarfgp
     public :: la_clarft
     public :: la_clarfx
     public :: la_clarfy
     public :: la_zlarf
     public :: la_zlarfb
     public :: la_zlarfg
     public :: la_zlarfgp
     public :: la_zlarft
     public :: la_zlarfx
     public :: la_zlarfy
#ifdef LA_WITH_XDP
     public :: la_ylarf
     public :: la_ylarfb
     public :: la_ylarfg
     public :: la_ylarfgp
     public :: la_ylarft
     public :: la_ylarfx
     public :: la_ylarfy
#endif
#ifdef LA_WITH_QP
     public :: la_wlarf
     public :: la_wlarfb
     public :: la_wlarfg
     public :: la_wlarfgp
     public :: la_wlarft
     public :: la_wlarfx
     public :: la_wlarfy
#endif

     contains

     !> SLARF: applies a real elementary reflector H to a real m by n matrix
     !> C, from either the left or the right. H is represented in the form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.

     pure subroutine la_slarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           real(sp),intent(in) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(in) :: v(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= zero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-zero row in v.
              do while (lastv > 0 .and. v(i) == zero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-zero column in c(1:lastv,:).
                 lastc = la_ilaslc(lastv,n,c,ldc)
              else
           ! scan for the last non-zero row in c(:,1:lastv).
                 lastc = la_ilaslr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_sp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
                 call la_sgemv('TRANSPOSE',lastv,lastc,one,c,ldc,v,incv,zero,work,1 &
                           )
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
                 call la_sger(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_sgemv('NO TRANSPOSE',lastc,lastv,one,c,ldc,v,incv,zero,work, &
                            1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
                 call la_sger(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_slarf
     !> DLARF: applies a real elementary reflector H to a real m by n matrix
     !> C, from either the left or the right. H is represented in the form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.

     pure subroutine la_dlarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           real(dp),intent(in) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(in) :: v(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= zero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-zero row in v.
              do while (lastv > 0 .and. v(i) == zero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-zero column in c(1:lastv,:).
                 lastc = la_iladlc(lastv,n,c,ldc)
              else
           ! scan for the last non-zero row in c(:,1:lastv).
                 lastc = la_iladlr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_dp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
                 call la_dgemv('TRANSPOSE',lastv,lastc,one,c,ldc,v,incv,zero,work,1 &
                           )
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
                 call la_dger(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_dgemv('NO TRANSPOSE',lastc,lastv,one,c,ldc,v,incv,zero,work, &
                            1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
                 call la_dger(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_dlarf
#ifdef LA_WITH_XDP
     !> XLARF: applies a real elementary reflector H to a real m by n matrix
     !> C, from either the left or the right. H is represented in the form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.

     pure subroutine la_xlarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           real(xdp),intent(in) :: tau
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(in) :: v(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= zero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-zero row in v.
              do while (lastv > 0 .and. v(i) == zero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-zero column in c(1:lastv,:).
                 lastc = la_ilaxlc(lastv,n,c,ldc)
              else
           ! scan for the last non-zero row in c(:,1:lastv).
                 lastc = la_ilaxlr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_xdp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
                 call la_xgemv('TRANSPOSE',lastv,lastc,one,c,ldc,v,incv,zero,work,1 &
                           )
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
                 call la_xger(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_xgemv('NO TRANSPOSE',lastc,lastv,one,c,ldc,v,incv,zero,work, &
                            1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
                 call la_xger(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_xlarf
#endif
#ifdef LA_WITH_QP
     !> QLARF: applies a real elementary reflector H to a real m by n matrix
     !> C, from either the left or the right. H is represented in the form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix.

     pure subroutine la_qlarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           real(qp),intent(in) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(in) :: v(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= zero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-zero row in v.
              do while (lastv > 0 .and. v(i) == zero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-zero column in c(1:lastv,:).
                 lastc = la_ilaqlc(lastv,n,c,ldc)
              else
           ! scan for the last non-zero row in c(:,1:lastv).
                 lastc = la_ilaqlr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_qp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**t * v(1:lastv,1)
                 call la_qgemv('TRANSPOSE',lastv,lastc,one,c,ldc,v,incv,zero,work,1 &
                           )
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**t
                 call la_qger(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_qgemv('NO TRANSPOSE',lastc,lastv,one,c,ldc,v,incv,zero,work, &
                            1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**t
                 call la_qger(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_qlarf
#endif

     !> SLARFB: applies a real block reflector H or its transpose H**T to a
     !> real m by n matrix C, from either the left or the right.

     pure subroutine la_slarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(in) :: t(ldt,*),v(ldv,*)
           real(sp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_scopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_strmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2
                       call la_sgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                  ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_strmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c2 := c2 - v2 * w**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(k + 1,1) &
                                 ,ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_strmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_scopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_strmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_strmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c2 := c2 - w * v2**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(k + 1,1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_strmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_scopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_strmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v(m - k + &
                              1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1
                       call la_sgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_strmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c1 := c1 - v1 * w**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv, &
                                 work,ldwork,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_strmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v(m - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h'  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_scopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_strmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v(n - k + &
                              1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c,ldc, &
                                 v,ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_strmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c1 := c1 - w * v1**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_strmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v(n - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_scopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_strmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2**t
                       call la_sgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                 ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_strmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c2 := c2 - v2**t * w**t
                       call la_sgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(1,k + 1), &
                                 ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_strmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_scopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_strmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c(1,k + 1), &
                                  ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_strmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(1,k + 1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_strmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_scopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_strmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v(1,m - k + &
                              1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1**t
                       call la_sgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c,ldc,v,ldv, &
                                  one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_strmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c1 := c1 - v1**t * w**t
                       call la_sgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv,work, &
                                 ldwork,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_strmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_scopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_strmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v(1,n - k + &
                              1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**t
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_strmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_strmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_slarfb
     !> DLARFB: applies a real block reflector H or its transpose H**T to a
     !> real m by n matrix C, from either the left or the right.

     pure subroutine la_dlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(in) :: t(ldt,*),v(ldv,*)
           real(dp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_dcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_dtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2
                       call la_dgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                  ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_dtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c2 := c2 - v2 * w**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(k + 1,1) &
                                 ,ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_dtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_dcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_dtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_dtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c2 := c2 - w * v2**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(k + 1,1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_dtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_dcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_dtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v(m - k + &
                              1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1
                       call la_dgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_dtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c1 := c1 - v1 * w**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv, &
                                 work,ldwork,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_dtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v(m - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_dcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_dtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v(n - k + &
                              1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c,ldc, &
                                 v,ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_dtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c1 := c1 - w * v1**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_dtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v(n - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_dcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_dtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2**t
                       call la_dgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                 ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_dtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c2 := c2 - v2**t * w**t
                       call la_dgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(1,k + 1), &
                                 ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_dtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_dcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_dtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c(1,k + 1), &
                                  ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_dtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(1,k + 1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_dtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_dcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_dtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v(1,m - k + &
                              1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1**t
                       call la_dgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c,ldc,v,ldv, &
                                  one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_dtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c1 := c1 - v1**t * w**t
                       call la_dgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv,work, &
                                 ldwork,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_dtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h'  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_dcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_dtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v(1,n - k + &
                              1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**t
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_dtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_dtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dlarfb
#ifdef LA_WITH_XDP
     !> XLARFB: applies a real block reflector H or its transpose H**T to a
     !> real m by n matrix C, from either the left or the right.

     pure subroutine la_xlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(in) :: t(ldt,*),v(ldv,*)
           real(xdp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_xcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_xtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2
                       call la_xgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                  ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_xtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c2 := c2 - v2 * w**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(k + 1,1) &
                                 ,ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_xtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_xcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_xtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_xtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c2 := c2 - w * v2**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(k + 1,1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_xtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_xcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_xtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v(m - k + &
                              1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1
                       call la_xgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_xtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c1 := c1 - v1 * w**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv, &
                                 work,ldwork,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_xtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v(m - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_xcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_xtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v(n - k + &
                              1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c,ldc, &
                                 v,ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_xtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c1 := c1 - w * v1**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_xtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v(n - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_xcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_xtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2**t
                       call la_xgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                 ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_xtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c2 := c2 - v2**t * w**t
                       call la_xgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(1,k + 1), &
                                 ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_xtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_xcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_xtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c(1,k + 1), &
                                  ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_xtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(1,k + 1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_xtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_xcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_xtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v(1,m - k + &
                              1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1**t
                       call la_xgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c,ldc,v,ldv, &
                                  one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_xtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c1 := c1 - v1**t * w**t
                       call la_xgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv,work, &
                                 ldwork,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_xtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h'  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_xcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_xtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v(1,n - k + &
                              1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**t
                       call la_xgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_xtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_xtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_xlarfb
#endif
#ifdef LA_WITH_QP
     !> QLARFB: applies a real block reflector H or its transpose H**T to a
     !> real m by n matrix C, from either the left or the right.

     pure subroutine la_qlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(in) :: t(ldt,*),v(ldv,*)
           real(qp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'T'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_qcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_qtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2
                       call la_qgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                  ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_qtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c2 := c2 - v2 * w**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(k + 1,1) &
                                 ,ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_qtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_qcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_qtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_qtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c2 := c2 - w * v2**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(k + 1,1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**t
                    call la_qtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_qcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_qtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v(m - k + &
                              1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1
                       call la_qgemm('TRANSPOSE','NO TRANSPOSE',n,k,m - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_qtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**t
                    if (m > k) then
                       ! c1 := c1 - v1 * w**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv, &
                                 work,ldwork,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_qtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v(m - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_qcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_qtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v(n - k + &
                              1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,one,c,ldc, &
                                 v,ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_qtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**t
                    if (n > k) then
                       ! c1 := c1 - w * v1**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2**t
                    call la_qtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v(n - k + 1, &
                              1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c1**t
                    do j = 1,k
                       call la_qcopy(n,c(j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_qtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',n,k,one,v,ldv, &
                              work,ldwork)
                    if (m > k) then
                       ! w := w + c2**t * v2**t
                       call la_qgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c(k + 1,1), &
                                 ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_qtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c2 := c2 - v2**t * w**t
                       call la_qgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v(1,k + 1), &
                                 ldv,work,ldwork,one,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_qtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_qcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**t
                    call la_qtrmm('RIGHT','UPPER','TRANSPOSE','UNIT',m,k,one,v,ldv, &
                              work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c(1,k + 1), &
                                  ldc,v(1,k + 1),ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_qtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v(1,k + 1),ldv,one,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_qtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,one,v,ldv, &
                               work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**t * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    ! w := c2**t
                    do j = 1,k
                       call la_qcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_qtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',n,k,one,v(1,m - k + &
                              1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**t * v1**t
                       call la_qgemm('TRANSPOSE','TRANSPOSE',n,k,m - k,one,c,ldc,v,ldv, &
                                  one,work,ldwork)
                    end if
                    ! w := w * t**t  or  w * t
                    call la_qtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - v**t * w**t
                    if (m > k) then
                       ! c1 := c1 - v1**t * w**t
                       call la_qgemm('TRANSPOSE','TRANSPOSE',m - k,n,k,-one,v,ldv,work, &
                                 ldwork,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_qtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,one,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**t
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - work(i,j)
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h'  where  c = ( c1  c2 )
                    ! w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_qcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**t
                    call la_qtrmm('RIGHT','LOWER','TRANSPOSE','UNIT',m,k,one,v(1,n - k + &
                              1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**t
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',m,k,n - k,one,c,ldc,v, &
                                 ldv,one,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**t
                    call la_qtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,one,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-one,work, &
                                 ldwork,v,ldv,one,c,ldc)
                    end if
                    ! w := w * v2
                    call la_qtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,one,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qlarfb
#endif

     !> SLARFT: forms the triangular factor T of a real block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V

     pure subroutine la_slarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(sp),intent(out) :: t(ldt,*)
           real(sp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(i,prevlastv)
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(i,j)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**t * v(i:j,i)
                       call la_sgemv('TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1),ldv,v(i + &
                                 1,i),1,one,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**t
                       call la_sgemv('NO TRANSPOSE',i - 1,j - i,-tau(i),v(1,i + 1),ldv,v( &
                                  i,i + 1),ldv,one,t(1,i),1)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_strmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                 end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(n - k + i,j)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**t * v(j:n-k+i,i)
                          call la_sgemv('TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j,i + 1), &
                                    ldv,v(j,i),1,one,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**t
                          call la_sgemv('NO TRANSPOSE',k - i,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,one,t(i + 1,i),1)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_strmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_slarft
     !> DLARFT: forms the triangular factor T of a real block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V

     pure subroutine la_dlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(dp),intent(out) :: t(ldt,*)
           real(dp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(i,prevlastv)
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(i,j)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**t * v(i:j,i)
                       call la_dgemv('TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1),ldv,v(i + &
                                 1,i),1,one,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**t
                       call la_dgemv('NO TRANSPOSE',i - 1,j - i,-tau(i),v(1,i + 1),ldv,v( &
                                  i,i + 1),ldv,one,t(1,i),1)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_dtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                 end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(n - k + i,j)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**t * v(j:n-k+i,i)
                          call la_dgemv('TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j,i + 1), &
                                    ldv,v(j,i),1,one,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**t
                          call la_dgemv('NO TRANSPOSE',k - i,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,one,t(i + 1,i),1)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_dtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_dlarft
#ifdef LA_WITH_XDP
     !> XLARFT: forms the triangular factor T of a real block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V

     pure subroutine la_xlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(xdp),intent(out) :: t(ldt,*)
           real(xdp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(i,prevlastv)
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(i,j)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**t * v(i:j,i)
                       call la_xgemv('TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1),ldv,v(i + &
                                 1,i),1,one,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**t
                       call la_xgemv('NO TRANSPOSE',i - 1,j - i,-tau(i),v(1,i + 1),ldv,v( &
                                  i,i + 1),ldv,one,t(1,i),1)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_xtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                 end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(n - k + i,j)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**t * v(j:n-k+i,i)
                          call la_xgemv('TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j,i + 1), &
                                    ldv,v(j,i),1,one,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**t
                          call la_xgemv('NO TRANSPOSE',k - i,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,one,t(i + 1,i),1)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_xtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_xlarft
#endif
#ifdef LA_WITH_QP
     !> QLARFT: forms the triangular factor T of a real block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**T
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**T * T * V

     pure subroutine la_qlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           real(qp),intent(out) :: t(ldt,*)
           real(qp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(i,prevlastv)
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(i,j)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**t * v(i:j,i)
                       call la_qgemv('TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1),ldv,v(i + &
                                 1,i),1,one,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= zero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**t
                       call la_qgemv('NO TRANSPOSE',i - 1,j - i,-tau(i),v(1,i + 1),ldv,v( &
                                  i,i + 1),ldv,one,t(1,i),1)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_qtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                 end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == zero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = zero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(n - k + i,j)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**t * v(j:n-k+i,i)
                          call la_qgemv('TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j,i + 1), &
                                    ldv,v(j,i),1,one,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= zero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**t
                          call la_qgemv('NO TRANSPOSE',k - i,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,one,t(i + 1,i),1)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_qtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_qlarft
#endif

     !> SLARFX: applies a real elementary reflector H to a real m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_slarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           real(sp),intent(in) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(in) :: v(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           real(sp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5,v6, &
                     v7,v8,v9
           ! Executable Statements
           if (tau == zero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_slarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_slarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
410 return
     end subroutine la_slarfx
     !> DLARFX: applies a real elementary reflector H to a real m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_dlarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           real(dp),intent(in) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(in) :: v(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           real(dp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5,v6, &
                     v7,v8,v9
           ! Executable Statements
           if (tau == zero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_dlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_dlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_dlarfx
#ifdef LA_WITH_XDP
     !> XLARFX: applies a real elementary reflector H to a real m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_xlarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           real(xdp),intent(in) :: tau
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(in) :: v(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           real(xdp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5,v6, &
                     v7,v8,v9
           ! Executable Statements
           if (tau == zero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_xlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_xlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_xlarfx
#endif
#ifdef LA_WITH_QP
     !> QLARFX: applies a real elementary reflector H to a real m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**T
     !> where tau is a real scalar and v is a real vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_qlarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           real(qp),intent(in) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(in) :: v(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           real(qp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5,v6, &
                     v7,v8,v9
           ! Executable Statements
           if (tau == zero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_qlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_qlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = one - tau*v(1)*v(1)
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*v1
              v2 = v(2)
              t2 = tau*v2
              v3 = v(3)
              t3 = tau*v3
              v4 = v(4)
              t4 = tau*v4
              v5 = v(5)
              t5 = tau*v5
              v6 = v(6)
              t6 = tau*v6
              v7 = v(7)
              t7 = tau*v7
              v8 = v(8)
              t8 = tau*v8
              v9 = v(9)
              t9 = tau*v9
              v10 = v(10)
              t10 = tau*v10
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_qlarfx
#endif

     !> SLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n symmetric matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_slarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           real(sp),intent(in) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*)
           real(sp),intent(in) :: v(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: alpha
           ! Executable Statements
           if (tau == zero) return
           ! form  w:= c * v
           call la_ssymv(uplo,n,one,c,ldc,v,incv,zero,work,1)
           alpha = -half*tau*la_sdot(n,work,1,v,incv)
           call la_saxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_ssyr2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_slarfy
     !> DLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n symmetric matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_dlarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           real(dp),intent(in) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*)
           real(dp),intent(in) :: v(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: alpha
           ! Executable Statements
           if (tau == zero) return
           ! form  w:= c * v
           call la_dsymv(uplo,n,one,c,ldc,v,incv,zero,work,1)
           alpha = -half*tau*la_ddot(n,work,1,v,incv)
           call la_daxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_dsyr2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_dlarfy
#ifdef LA_WITH_XDP
     !> XLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n symmetric matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_xlarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           real(xdp),intent(in) :: tau
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*)
           real(xdp),intent(in) :: v(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: alpha
           ! Executable Statements
           if (tau == zero) return
           ! form  w:= c * v
           call la_xsymv(uplo,n,one,c,ldc,v,incv,zero,work,1)
           alpha = -half*tau*la_xdot(n,work,1,v,incv)
           call la_xaxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_xsyr2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_xlarfy
#endif
#ifdef LA_WITH_QP
     !> QLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n symmetric matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_qlarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           real(qp),intent(in) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*)
           real(qp),intent(in) :: v(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: alpha
           ! Executable Statements
           if (tau == zero) return
           ! form  w:= c * v
           call la_qsymv(uplo,n,one,c,ldc,v,incv,zero,work,1)
           alpha = -half*tau*la_qdot(n,work,1,v,incv)
           call la_qaxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_qsyr2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_qlarfy
#endif

     !> SLARFG: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, and x is an (n-1)-element real
     !> vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.
     !> Otherwise  1 <= tau <= 2.

     pure subroutine la_slarfg(n,alpha,x,incx,tau)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(sp),intent(inout) :: alpha
           real(sp),intent(out) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(sp) :: beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 1) then
              tau = zero
              return
           end if
           xnorm = la_snrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_slapy2(alpha,xnorm),alpha)
              safmin = la_slamch('S')/la_slamch('E')
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 rsafmn = one/safmin
                 10 continue
                 knt = knt + 1
                 call la_sscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alpha = alpha*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_snrm2(n - 1,x,incx)
                 beta = -sign(la_slapy2(alpha,xnorm),alpha)
              end if
              tau = (beta - alpha)/beta
              call la_sscal(n - 1,one/(alpha - beta),x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_slarfg
     !> DLARFG: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, and x is an (n-1)-element real
     !> vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.
     !> Otherwise  1 <= tau <= 2.

     pure subroutine la_dlarfg(n,alpha,x,incx,tau)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(dp),intent(inout) :: alpha
           real(dp),intent(out) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(dp) :: beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 1) then
              tau = zero
              return
           end if
           xnorm = la_dnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_dlapy2(alpha,xnorm),alpha)
              safmin = la_dlamch('S')/la_dlamch('E')
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 rsafmn = one/safmin
                 10 continue
                 knt = knt + 1
                 call la_dscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alpha = alpha*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_dnrm2(n - 1,x,incx)
                 beta = -sign(la_dlapy2(alpha,xnorm),alpha)
              end if
              tau = (beta - alpha)/beta
              call la_dscal(n - 1,one/(alpha - beta),x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_dlarfg
#ifdef LA_WITH_XDP
     !> XLARFG: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, and x is an (n-1)-element real
     !> vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.
     !> Otherwise  1 <= tau <= 2.

     pure subroutine la_xlarfg(n,alpha,x,incx,tau)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(xdp),intent(inout) :: alpha
           real(xdp),intent(out) :: tau
           ! Array Arguments
           real(xdp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(xdp) :: beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 1) then
              tau = zero
              return
           end if
           xnorm = la_xnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_xlapy2(alpha,xnorm),alpha)
              safmin = la_xlamch('S')/la_xlamch('E')
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 rsafmn = one/safmin
                 10 continue
                 knt = knt + 1
                 call la_xscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alpha = alpha*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_xnrm2(n - 1,x,incx)
                 beta = -sign(la_xlapy2(alpha,xnorm),alpha)
              end if
              tau = (beta - alpha)/beta
              call la_xscal(n - 1,one/(alpha - beta),x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_xlarfg
#endif
#ifdef LA_WITH_QP
     !> QLARFG: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, and x is an (n-1)-element real
     !> vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.
     !> Otherwise  1 <= tau <= 2.

     pure subroutine la_qlarfg(n,alpha,x,incx,tau)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(qp),intent(inout) :: alpha
           real(qp),intent(out) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(qp) :: beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 1) then
              tau = zero
              return
           end if
           xnorm = la_qnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_qlapy2(alpha,xnorm),alpha)
              safmin = la_qlamch('S')/la_qlamch('E')
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 rsafmn = one/safmin
                 10 continue
                 knt = knt + 1
                 call la_qscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alpha = alpha*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_qnrm2(n - 1,x,incx)
                 beta = -sign(la_qlapy2(alpha,xnorm),alpha)
              end if
              tau = (beta - alpha)/beta
              call la_qscal(n - 1,one/(alpha - beta),x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_qlarfg
#endif

     !> SLARFGP: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is non-negative, and x is
     !> an (n-1)-element real vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.

     subroutine la_slarfgp(n,alpha,x,incx,tau)
        use la_constants_sp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(sp),intent(inout) :: alpha
           real(sp),intent(out) :: tau
           ! Array Arguments
           real(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(sp) :: beta,bignum,savealpha,smlnum,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_snrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  [+/-1, 0; i], sign chosen so alpha >= 0.
              if (alpha >= zero) then
                 ! when tau.eq.zero, the vector is special-cased to be
                 ! all zeros in the application routines.  we do not need
                 ! to clear it.
                 tau = zero
              else
                 ! however, the application routines rely on explicit
                 ! zero checks when tau.ne.zero, and we must clear x.
                 tau = two
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = 0
                 end do
                 alpha = -alpha
              end if
           else
              ! general case
              beta = sign(la_slapy2(alpha,xnorm),alpha)
              smlnum = la_slamch('S')/la_slamch('E')
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 bignum = one/smlnum
                 10 continue
                 knt = knt + 1
                 call la_sscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alpha = alpha*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_snrm2(n - 1,x,incx)
                 beta = sign(la_slapy2(alpha,xnorm),alpha)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alpha = xnorm*(xnorm/alpha)
                 tau = alpha/beta
                 alpha = -alpha
              end if
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero. this explains the next if statement.
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 if (savealpha >= zero) then
                    tau = zero
                 else
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = 0
                    end do
                    beta = -savealpha
                 end if
              else
                 ! this is the general case.
                 call la_sscal(n - 1,one/alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_slarfgp
     !> DLARFGP: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is non-negative, and x is
     !> an (n-1)-element real vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.

     subroutine la_dlarfgp(n,alpha,x,incx,tau)
        use la_constants_dp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(dp),intent(inout) :: alpha
           real(dp),intent(out) :: tau
           ! Array Arguments
           real(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(dp) :: beta,bignum,savealpha,smlnum,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_dnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  [+/-1, 0; i], sign chosen so alpha >= 0
              if (alpha >= zero) then
                 ! when tau.eq.zero, the vector is special-cased to be
                 ! all zeros in the application routines.  we do not need
                 ! to clear it.
                 tau = zero
              else
                 ! however, the application routines rely on explicit
                 ! zero checks when tau.ne.zero, and we must clear x.
                 tau = two
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = 0
                 end do
                 alpha = -alpha
              end if
           else
              ! general case
              beta = sign(la_dlapy2(alpha,xnorm),alpha)
              smlnum = la_dlamch('S')/la_dlamch('E')
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 bignum = one/smlnum
                 10 continue
                 knt = knt + 1
                 call la_dscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alpha = alpha*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_dnrm2(n - 1,x,incx)
                 beta = sign(la_dlapy2(alpha,xnorm),alpha)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alpha = xnorm*(xnorm/alpha)
                 tau = alpha/beta
                 alpha = -alpha
              end if
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero. this explains the next if statement.
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 if (savealpha >= zero) then
                    tau = zero
                 else
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = 0
                    end do
                    beta = -savealpha
                 end if
              else
                 ! this is the general case.
                 call la_dscal(n - 1,one/alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_dlarfgp
#ifdef LA_WITH_XDP
     !> XLARFGP: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is non-negative, and x is
     !> an (n-1)-element real vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.

     subroutine la_xlarfgp(n,alpha,x,incx,tau)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(xdp),intent(inout) :: alpha
           real(xdp),intent(out) :: tau
           ! Array Arguments
           real(xdp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(xdp) :: beta,bignum,savealpha,smlnum,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_xnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  [+/-1, 0; i], sign chosen so alpha >= 0
              if (alpha >= zero) then
                 ! when tau.eq.zero, the vector is special-cased to be
                 ! all zeros in the application routines.  we do not need
                 ! to clear it.
                 tau = zero
              else
                 ! however, the application routines rely on explicit
                 ! zero checks when tau.ne.zero, and we must clear x.
                 tau = two
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = 0
                 end do
                 alpha = -alpha
              end if
           else
              ! general case
              beta = sign(la_xlapy2(alpha,xnorm),alpha)
              smlnum = la_xlamch('S')/la_xlamch('E')
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 bignum = one/smlnum
                 10 continue
                 knt = knt + 1
                 call la_xscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alpha = alpha*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_xnrm2(n - 1,x,incx)
                 beta = sign(la_xlapy2(alpha,xnorm),alpha)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alpha = xnorm*(xnorm/alpha)
                 tau = alpha/beta
                 alpha = -alpha
              end if
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero. this explains the next if statement.
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 if (savealpha >= zero) then
                    tau = zero
                 else
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = 0
                    end do
                    beta = -savealpha
                 end if
              else
                 ! this is the general case.
                 call la_xscal(n - 1,one/alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_xlarfgp
#endif
#ifdef LA_WITH_QP
     !> QLARFGP: generates a real elementary reflector H of order n, such
     !> that
     !> H * ( alpha ) = ( beta ),   H**T * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is non-negative, and x is
     !> an (n-1)-element real vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**T ) ,
     !> ( v )
     !> where tau is a real scalar and v is a real (n-1)-element
     !> vector.
     !> If the elements of x are all zero, then tau = 0 and H is taken to be
     !> the unit matrix.

     subroutine la_qlarfgp(n,alpha,x,incx,tau)
        use la_constants_qp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(qp),intent(inout) :: alpha
           real(qp),intent(out) :: tau
           ! Array Arguments
           real(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(qp) :: beta,bignum,savealpha,smlnum,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_qnrm2(n - 1,x,incx)
           if (xnorm == zero) then
              ! h  =  [+/-1, 0; i], sign chosen so alpha >= 0
              if (alpha >= zero) then
                 ! when tau.eq.zero, the vector is special-cased to be
                 ! all zeros in the application routines.  we do not need
                 ! to clear it.
                 tau = zero
              else
                 ! however, the application routines rely on explicit
                 ! zero checks when tau.ne.zero, and we must clear x.
                 tau = two
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = 0
                 end do
                 alpha = -alpha
              end if
           else
              ! general case
              beta = sign(la_qlapy2(alpha,xnorm),alpha)
              smlnum = la_qlamch('S')/la_qlamch('E')
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 bignum = one/smlnum
                 10 continue
                 knt = knt + 1
                 call la_qscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alpha = alpha*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_qnrm2(n - 1,x,incx)
                 beta = sign(la_qlapy2(alpha,xnorm),alpha)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alpha = xnorm*(xnorm/alpha)
                 tau = alpha/beta
                 alpha = -alpha
              end if
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero. this explains the next if statement.
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 if (savealpha >= zero) then
                    tau = zero
                 else
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = 0
                    end do
                    beta = -savealpha
                 end if
              else
                 ! this is the general case.
                 call la_qscal(n - 1,one/alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_qlarfgp
#endif

     !> CLARF: applies a complex elementary reflector H to a complex M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H (the conjugate transpose of H), supply conjg(tau) instead
     !> tau.

     pure subroutine la_clarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           complex(sp),intent(in) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(in) :: v(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= czero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-czero row in v.
              do while (lastv > 0 .and. v(i) == czero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-czero column in c(1:lastv,:).
                 lastc = la_ilaclc(lastv,n,c,ldc)
              else
           ! scan for the last non-czero row in c(:,1:lastv).
                 lastc = la_ilaclr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_sp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**h * v(1:lastv,1)
                 call la_cgemv('CONJUGATE TRANSPOSE',lastv,lastc,cone,c,ldc,v,incv, &
                           czero,work,1)
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**h
                 call la_cgerc(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_cgemv('NO TRANSPOSE',lastc,lastv,cone,c,ldc,v,incv,czero, &
                           work,1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**h
                 call la_cgerc(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_clarf
     !> ZLARF: applies a complex elementary reflector H to a complex M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H, supply conjg(tau) instead
     !> tau.

     pure subroutine la_zlarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           complex(dp),intent(in) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(in) :: v(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= czero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-czero row in v.
              do while (lastv > 0 .and. v(i) == czero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-czero column in c(1:lastv,:).
                 lastc = la_ilazlc(lastv,n,c,ldc)
              else
           ! scan for the last non-czero row in c(:,1:lastv).
                 lastc = la_ilazlr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_dp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**h * v(1:lastv,1)
                 call la_zgemv('CONJUGATE TRANSPOSE',lastv,lastc,cone,c,ldc,v,incv, &
                           czero,work,1)
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**h
                 call la_zgerc(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_zgemv('NO TRANSPOSE',lastc,lastv,cone,c,ldc,v,incv,czero, &
                           work,1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**h
                 call la_zgerc(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_zlarf
#ifdef LA_WITH_XDP
     !> YLARF: applies a complex elementary reflector H to a complex M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H, supply conjg(tau) instead
     !> tau.

     pure subroutine la_ylarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           complex(xdp),intent(in) :: tau
           ! Array Arguments
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(in) :: v(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= czero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-czero row in v.
              do while (lastv > 0 .and. v(i) == czero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-czero column in c(1:lastv,:).
                 lastc = la_ilaylc(lastv,n,c,ldc)
              else
           ! scan for the last non-czero row in c(:,1:lastv).
                 lastc = la_ilaylr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_xdp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**h * v(1:lastv,1)
                 call la_ygemv('CONJUGATE TRANSPOSE',lastv,lastc,cone,c,ldc,v,incv, &
                           czero,work,1)
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**h
                 call la_ygerc(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_ygemv('NO TRANSPOSE',lastc,lastv,cone,c,ldc,v,incv,czero, &
                           work,1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**h
                 call la_ygerc(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_ylarf
#endif
#ifdef LA_WITH_QP
     !> WLARF: applies a complex elementary reflector H to a complex M-by-N
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix.
     !> To apply H**H, supply conjg(tau) instead
     !> tau.

     pure subroutine la_wlarf(side,m,n,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: incv,ldc,m,n
           complex(qp),intent(in) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(in) :: v(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: applyleft
           integer(ilp) :: i,lastv,lastc
           ! Executable Statements
           applyleft = la_lsame(side,'L')
           lastv = 0
           lastc = 0
           if (tau /= czero) then
           ! set up variables for scanning v.  lastv begins pointing to the end
           ! of v.
              if (applyleft) then
                 lastv = m
              else
                 lastv = n
              end if
              if (incv > 0) then
                 i = 1 + (lastv - 1)*incv
              else
                 i = 1
              end if
           ! look for the last non-czero row in v.
              do while (lastv > 0 .and. v(i) == czero)
                 lastv = lastv - 1
                 i = i - incv
              end do
              if (applyleft) then
           ! scan for the last non-czero column in c(1:lastv,:).
                 lastc = la_ilawlc(lastv,n,c,ldc)
              else
           ! scan for the last non-czero row in c(:,1:lastv).
                 lastc = la_ilawlr(m,lastv,c,ldc)
              end if
           end if
           ! note that lastc.eq.0_qp renders the blas operations null; no special
           ! case is needed at this level.
           if (applyleft) then
              ! form  h * c
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastv,1:lastc)**h * v(1:lastv,1)
                 call la_wgemv('CONJUGATE TRANSPOSE',lastv,lastc,cone,c,ldc,v,incv, &
                           czero,work,1)
                 ! c(1:lastv,1:lastc) := c(...) - v(1:lastv,1) * w(1:lastc,1)**h
                 call la_wgerc(lastv,lastc,-tau,v,incv,work,1,c,ldc)
              end if
           else
              ! form  c * h
              if (lastv > 0) then
                 ! w(1:lastc,1) := c(1:lastc,1:lastv) * v(1:lastv,1)
                 call la_wgemv('NO TRANSPOSE',lastc,lastv,cone,c,ldc,v,incv,czero, &
                           work,1)
                 ! c(1:lastc,1:lastv) := c(...) - w(1:lastc,1) * v(1:lastv,1)**h
                 call la_wgerc(lastc,lastv,-tau,work,1,v,incv,c,ldc)
              end if
           end if
           return
     end subroutine la_wlarf
#endif

     !> CLARFB: applies a complex block reflector H or its transpose H**H to a
     !> complex M-by-N matrix C, from either the left or the right.

     pure subroutine la_clarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(in) :: t(ldt,*),v(ldv,*)
           complex(sp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_ccopy(n,c(j,1),ldc,work(1,j),1)
                       call la_clacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h *v2
                       call la_cgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c(k + 1,1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ctrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c2 := c2 - v2 * w**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v(k + 1,1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ctrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_ccopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ctrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c2 := c2 - w * v2**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v(k + 1,1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ctrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_ccopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_clacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v(m - &
                              k + 1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1
                       call la_cgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ctrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c1 := c1 - v1 * w**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(m - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_ccopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v(n - &
                              k + 1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c,ldc, &
                                 v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ctrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c1 := c1 - w * v1**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(n - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_ccopy(n,c(j,1),ldc,work(1,j),1)
                       call la_clacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2**h
                       call la_cgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c(k + 1,1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ctrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c2 := c2 - v2**h * w**h
                       call la_cgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v(1,k + 1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_ccopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c(1,k + 1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ctrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v(1,k + 1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_ctrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_ccopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_clacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ctrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(1,m - k + 1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1**h
                       call la_cgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ctrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c1 := c1 - v1**h * w**h
                       call la_cgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_ccopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ctrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(1,n - k + 1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**h
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ctrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ctrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_clarfb
     !> ZLARFB: applies a complex block reflector H or its transpose H**H to a
     !> complex M-by-N matrix C, from either the left or the right.

     pure subroutine la_zlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(in) :: t(ldt,*),v(ldv,*)
           complex(dp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_zcopy(n,c(j,1),ldc,work(1,j),1)
                       call la_zlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2
                       call la_zgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c(k + 1,1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ztrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c2 := c2 - v2 * w**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v(k + 1,1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ztrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_zcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ztrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c2 := c2 - w * v2**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v(k + 1,1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ztrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_zcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_zlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v(m - &
                              k + 1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1
                       call la_zgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ztrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c1 := c1 - v1 * w**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(m - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_zcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v(n - &
                              k + 1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c,ldc, &
                                 v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ztrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c1 := c1 - w * v1**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(n - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_zcopy(n,c(j,1),ldc,work(1,j),1)
                       call la_zlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2**h
                       call la_zgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c(k + 1,1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ztrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c2 := c2 - v2**h * w**h
                       call la_zgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v(1,k + 1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_zcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c(1,k + 1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ztrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v(1,k + 1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_ztrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_zcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_zlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ztrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(1,m - k + 1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1**h
                       call la_zgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ztrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c1 := c1 - v1**h * w**h
                       call la_zgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_zcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ztrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(1,n - k + 1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**h
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ztrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ztrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_zlarfb
#ifdef LA_WITH_XDP
     !> YLARFB: applies a complex block reflector H or its transpose H**H to a
     !> complex M-by-N matrix C, from either the left or the right.

     pure subroutine la_ylarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(in) :: t(ldt,*),v(ldv,*)
           complex(xdp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_ycopy(n,c(j,1),ldc,work(1,j),1)
                       call la_ylacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2
                       call la_ygemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c(k + 1,1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ytrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c2 := c2 - v2 * w**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v(k + 1,1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ytrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_ycopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ytrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c2 := c2 - w * v2**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v(k + 1,1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_ytrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_ycopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_ylacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v(m - &
                              k + 1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1
                       call la_ygemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ytrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c1 := c1 - v1 * w**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ytrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(m - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_ycopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v(n - &
                              k + 1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c,ldc, &
                                 v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ytrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c1 := c1 - w * v1**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_ytrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(n - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_ycopy(n,c(j,1),ldc,work(1,j),1)
                       call la_ylacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ytrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2**h
                       call la_ygemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c(k + 1,1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ytrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c2 := c2 - v2**h * w**h
                       call la_ygemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v(1,k + 1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_ycopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_ytrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c(1,k + 1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ytrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v(1,k + 1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_ytrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_ycopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_ylacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ytrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(1,m - k + 1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1**h
                       call la_ygemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_ytrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c1 := c1 - v1**h * w**h
                       call la_ygemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_ycopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_ytrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(1,n - k + 1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**h
                       call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_ytrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_ytrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ylarfb
#endif
#ifdef LA_WITH_QP
     !> WLARFB: applies a complex block reflector H or its transpose H**H to a
     !> complex M-by-N matrix C, from either the left or the right.

     pure subroutine la_wlarfb(side,trans,direct,storev,m,n,k,v,ldv,t,ldt,c,ldc, &
               work,ldwork)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,side,storev,trans
           integer(ilp),intent(in) :: k,ldc,ldt,ldv,ldwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(in) :: t(ldt,*),v(ldv,*)
           complex(qp),intent(out) :: work(ldwork,*)
        ! =====================================================================

           ! Local Scalars
           character :: transt
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ! quick return if possible
           if (m <= 0 .or. n <= 0) return
           if (la_lsame(trans,'N')) then
              transt = 'C'
           else
              transt = 'N'
           end if
           if (la_lsame(storev,'C')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1 )    (first k rows)
                           ! ( v2 )
                 ! where  v1  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_wcopy(n,c(j,1),ldc,work(1,j),1)
                       call la_wlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2
                       call la_wgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c(k + 1,1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_wtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c2 := c2 - v2 * w**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v(k + 1,1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_wtrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_wcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1
                    call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2
                       call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c(1,k + &
                                 1),ldc,v(k + 1,1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_wtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c2 := c2 - w * v2**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v(k + 1,1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1**h
                    call la_wtrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1 )
                           ! ( v2 )    (last k rows)
                 ! where  v2  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v  =  (c1**h * v1 + c2**h * v2)  (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_wcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_wlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v(m - &
                              k + 1,1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1
                       call la_wgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',n,k,m - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_wtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v * w**h
                    if (m > k) then
                       ! c1 := c1 - v1 * w**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k,-cone, &
                                 v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(m - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_wcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2
                    call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v(n - &
                              k + 1,1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1
                       call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m,k,n - k,cone,c,ldc, &
                                 v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_wtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v**h
                    if (n > k) then
                       ! c1 := c1 - w * v1**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,n - k,k,-cone, &
                                 work,ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2**h
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(n - k + 1,1),ldv,work,ldwork)
                    ! c2 := c2 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           else if (la_lsame(storev,'R')) then
              if (la_lsame(direct,'F')) then
                 ! let  v =  ( v1  v2 )    (v1: first k columns)
                 ! where  v1  is unit upper triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c1**h
                    do j = 1,k
                       call la_wcopy(n,c(j,1),ldc,work(1,j),1)
                       call la_wlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v,ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c2**h * v2**h
                       call la_wgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c(k + 1,1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_wtrmm('RIGHT','UPPER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c2 := c2 - v2**h * w**h
                       call la_wgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v(1,k + 1),ldv,work,ldwork,cone,c(k + 1,1),ldc)
                    end if
                    ! w := w * v1
                    call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',n,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(j,i) = c(j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c1
                    do j = 1,k
                       call la_wcopy(m,c(1,j),1,work(1,j),1)
                    end do
                    ! w := w * v1**h
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v,ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c2 * v2**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c(1,k + 1),ldc,v(1,k + 1),ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_wtrmm('RIGHT','UPPER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c2 := c2 - w * v2
                       call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v(1,k + 1),ldv,cone,c(1,k + 1),ldc)
                    end if
                    ! w := w * v1
                    call la_wtrmm('RIGHT','UPPER','NO TRANSPOSE','UNIT',m,k,cone,v, &
                              ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,j) = c(i,j) - work(i,j)
                       end do
                    end do
                 end if
              else
                 ! let  v =  ( v1  v2 )    (v2: last k columns)
                 ! where  v2  is unit lower triangular.
                 if (la_lsame(side,'L')) then
                    ! form  h * c  or  h**h * c  where  c = ( c1 )
                                                          ! ( c2 )
                    ! w := c**h * v**h  =  (c1**h * v1**h + c2**h * v2**h) (stored in work)
                    ! w := c2**h
                    do j = 1,k
                       call la_wcopy(n,c(m - k + j,1),ldc,work(1,j),1)
                       call la_wlacgv(n,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_wtrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',n,k,cone, &
                               v(1,m - k + 1),ldv,work,ldwork)
                    if (m > k) then
                       ! w := w + c1**h * v1**h
                       call la_wgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',n,k,m - k, &
                                 cone,c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t**h  or  w * t
                    call la_wtrmm('RIGHT','LOWER',transt,'NON-UNIT',n,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - v**h * w**h
                    if (m > k) then
                       ! c1 := c1 - v1**h * w**h
                       call la_wgemm('CONJUGATE TRANSPOSE','CONJUGATE TRANSPOSE',m - k,n,k, &
                                 -cone,v,ldv,work,ldwork,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',n,k,cone,v(1, &
                              m - k + 1),ldv,work,ldwork)
                    ! c2 := c2 - w**h
                    do j = 1,k
                       do i = 1,n
                          c(m - k + j,i) = c(m - k + j,i) - conjg(work(i,j))
                       end do
                    end do
                 else if (la_lsame(side,'R')) then
                    ! form  c * h  or  c * h**h  where  c = ( c1  c2 )
                    ! w := c * v**h  =  (c1*v1**h + c2*v2**h)  (stored in work)
                    ! w := c2
                    do j = 1,k
                       call la_wcopy(m,c(1,n - k + j),1,work(1,j),1)
                    end do
                    ! w := w * v2**h
                    call la_wtrmm('RIGHT','LOWER','CONJUGATE TRANSPOSE','UNIT',m,k,cone, &
                               v(1,n - k + 1),ldv,work,ldwork)
                    if (n > k) then
                       ! w := w + c1 * v1**h
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m,k,n - k,cone, &
                                 c,ldc,v,ldv,cone,work,ldwork)
                    end if
                    ! w := w * t  or  w * t**h
                    call la_wtrmm('RIGHT','LOWER',trans,'NON-UNIT',m,k,cone,t,ldt, &
                              work,ldwork)
                    ! c := c - w * v
                    if (n > k) then
                       ! c1 := c1 - w * v1
                       call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m,n - k,k,-cone,work, &
                                 ldwork,v,ldv,cone,c,ldc)
                    end if
                    ! w := w * v2
                    call la_wtrmm('RIGHT','LOWER','NO TRANSPOSE','UNIT',m,k,cone,v(1, &
                              n - k + 1),ldv,work,ldwork)
                    ! c1 := c1 - w
                    do j = 1,k
                       do i = 1,m
                          c(i,n - k + j) = c(i,n - k + j) - work(i,j)
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wlarfb
#endif

     !> CLARFG: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, with beta real, and x is an
     !> (n-1)-element complex vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.
     !> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .

     pure subroutine la_clarfg(n,alpha,x,incx,tau)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(sp),intent(inout) :: alpha
           complex(sp),intent(out) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(sp) :: alphi,alphr,beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,real,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_scnrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=sp)
           alphi = aimag(alpha)
           if (xnorm == zero .and. alphi == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_slapy3(alphr,alphi,xnorm),alphr)
              safmin = la_slamch('S')/la_slamch('E')
              rsafmn = one/safmin
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_csscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alphi = alphi*rsafmn
                 alphr = alphr*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_scnrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=sp)
                 beta = -sign(la_slapy3(alphr,alphi,xnorm),alphr)
              end if
              tau = cmplx((beta - alphr)/beta,-alphi/beta,KIND=sp)
              alpha = la_cladiv(cmplx(one,KIND=sp),alpha - beta)
              call la_cscal(n - 1,alpha,x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_clarfg
     !> ZLARFG: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, with beta real, and x is an
     !> (n-1)-element complex vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.
     !> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .

     pure subroutine la_zlarfg(n,alpha,x,incx,tau)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(dp),intent(inout) :: alpha
           complex(dp),intent(out) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(dp) :: alphi,alphr,beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_dznrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=dp)
           alphi = aimag(alpha)
           if (xnorm == zero .and. alphi == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_dlapy3(alphr,alphi,xnorm),alphr)
              safmin = la_dlamch('S')/la_dlamch('E')
              rsafmn = one/safmin
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_zdscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alphi = alphi*rsafmn
                 alphr = alphr*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_dznrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=dp)
                 beta = -sign(la_dlapy3(alphr,alphi,xnorm),alphr)
              end if
              tau = cmplx((beta - alphr)/beta,-alphi/beta,KIND=dp)
              alpha = la_zladiv(cmplx(one,KIND=dp),alpha - beta)
              call la_zscal(n - 1,alpha,x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_zlarfg
#ifdef LA_WITH_XDP
     !> YLARFG: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, with beta real, and x is an
     !> (n-1)-element complex vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.
     !> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .

     pure subroutine la_ylarfg(n,alpha,x,incx,tau)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(xdp),intent(inout) :: alpha
           complex(xdp),intent(out) :: tau
           ! Array Arguments
           complex(xdp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(xdp) :: alphi,alphr,beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_xynrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=xdp)
           alphi = aimag(alpha)
           if (xnorm == zero .and. alphi == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_xlapy3(alphr,alphi,xnorm),alphr)
              safmin = la_xlamch('S')/la_xlamch('E')
              rsafmn = one/safmin
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_yxscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alphi = alphi*rsafmn
                 alphr = alphr*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_xynrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=xdp)
                 beta = -sign(la_xlapy3(alphr,alphi,xnorm),alphr)
              end if
              tau = cmplx((beta - alphr)/beta,-alphi/beta,KIND=xdp)
              alpha = la_yladiv(cmplx(one,KIND=xdp),alpha - beta)
              call la_yscal(n - 1,alpha,x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_ylarfg
#endif
#ifdef LA_WITH_QP
     !> WLARFG: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, with beta real, and x is an
     !> (n-1)-element complex vector. H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.
     !> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .

     pure subroutine la_wlarfg(n,alpha,x,incx,tau)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(qp),intent(inout) :: alpha
           complex(qp),intent(out) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(qp) :: alphi,alphr,beta,rsafmn,safmin,xnorm
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_qwnrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=qp)
           alphi = aimag(alpha)
           if (xnorm == zero .and. alphi == zero) then
              ! h  =  i
              tau = zero
           else
              ! general case
              beta = -sign(la_qlapy3(alphr,alphi,xnorm),alphr)
              safmin = la_qlamch('S')/la_qlamch('E')
              rsafmn = one/safmin
              knt = 0
              if (abs(beta) < safmin) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_wqscal(n - 1,rsafmn,x,incx)
                 beta = beta*rsafmn
                 alphi = alphi*rsafmn
                 alphr = alphr*rsafmn
                 if ((abs(beta) < safmin) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least safmin
                 xnorm = la_qwnrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=qp)
                 beta = -sign(la_qlapy3(alphr,alphi,xnorm),alphr)
              end if
              tau = cmplx((beta - alphr)/beta,-alphi/beta,KIND=qp)
              alpha = la_wladiv(cmplx(one,KIND=qp),alpha - beta)
              call la_wscal(n - 1,alpha,x,incx)
              ! if alpha is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*safmin
              end do
              alpha = beta
           end if
           return
     end subroutine la_wlarfg
#endif

     !> CLARFGP: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is real and non-negative, and
     !> x is an (n-1)-element complex vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.

     subroutine la_clarfgp(n,alpha,x,incx,tau)
        use la_constants_sp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(sp),intent(inout) :: alpha
           complex(sp),intent(out) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(sp) :: alphi,alphr,beta,bignum,smlnum,xnorm
           complex(sp) :: savealpha
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,real,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_scnrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=sp)
           alphi = aimag(alpha)
           if (xnorm == zero) then
              ! h  =  [1-alpha/abs(alpha) 0; 0 i], sign chosen so alpha >= 0.
              if (alphi == zero) then
                 if (alphr >= zero) then
                    ! when tau.eq.zero, the vector is special-cased to be
                    ! all zeros in the application routines.  we do not need
                    ! to clear it.
                    tau = zero
                 else
                    ! however, the application routines rely on explicit
                    ! zero checks when tau.ne.zero, and we must clear x.
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    alpha = -alpha
                 end if
              else
                 ! only "reflecting" the diagonal entry to be real and non-negative.
                 xnorm = la_slapy2(alphr,alphi)
                 tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=sp)
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = zero
                 end do
                 alpha = xnorm
              end if
           else
              ! general case
              beta = sign(la_slapy3(alphr,alphi,xnorm),alphr)
              smlnum = la_slamch('S')/la_slamch('E')
              bignum = one/smlnum
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_csscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alphi = alphi*bignum
                 alphr = alphr*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_scnrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=sp)
                 beta = sign(la_slapy3(alphr,alphi,xnorm),alphr)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alphr = alphi*(alphi/real(alpha,KIND=sp))
                 alphr = alphr + xnorm*(xnorm/real(alpha,KIND=sp))
                 tau = cmplx(alphr/beta,-alphi/beta,KIND=sp)
                 alpha = cmplx(-alphr,alphi,KIND=sp)
              end if
              alpha = la_cladiv(cmplx(one,KIND=sp),alpha)
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero (or two or whatever makes a nonnegative real number for beta).
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 alphr = real(savealpha,KIND=sp)
                 alphi = aimag(savealpha)
                 if (alphi == zero) then
                    if (alphr >= zero) then
                       tau = zero
                    else
                       tau = two
                       do j = 1,n - 1
                          x(1 + (j - 1)*incx) = zero
                       end do
                       beta = real(-savealpha,KIND=sp)
                    end if
                 else
                    xnorm = la_slapy2(alphr,alphi)
                    tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=sp)
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    beta = xnorm
                 end if
              else
                 ! this is the general case.
                 call la_cscal(n - 1,alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_clarfgp
     !> ZLARFGP: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is real and non-negative, and
     !> x is an (n-1)-element complex vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.

     subroutine la_zlarfgp(n,alpha,x,incx,tau)
        use la_constants_dp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(dp),intent(inout) :: alpha
           complex(dp),intent(out) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(dp) :: alphi,alphr,beta,bignum,smlnum,xnorm
           complex(dp) :: savealpha
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_dznrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=dp)
           alphi = aimag(alpha)
           if (xnorm == zero) then
              ! h  =  [1-alpha/abs(alpha) 0; 0 i], sign chosen so alpha >= 0.
              if (alphi == zero) then
                 if (alphr >= zero) then
                    ! when tau.eq.zero, the vector is special-cased to be
                    ! all zeros in the application routines.  we do not need
                    ! to clear it.
                    tau = zero
                 else
                    ! however, the application routines rely on explicit
                    ! zero checks when tau.ne.zero, and we must clear x.
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    alpha = -alpha
                 end if
              else
                 ! only "reflecting" the diagonal entry to be real and non-negative.
                 xnorm = la_dlapy2(alphr,alphi)
                 tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=dp)
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = zero
                 end do
                 alpha = xnorm
              end if
           else
              ! general case
              beta = sign(la_dlapy3(alphr,alphi,xnorm),alphr)
              smlnum = la_dlamch('S')/la_dlamch('E')
              bignum = one/smlnum
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_zdscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alphi = alphi*bignum
                 alphr = alphr*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_dznrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=dp)
                 beta = sign(la_dlapy3(alphr,alphi,xnorm),alphr)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alphr = alphi*(alphi/real(alpha,KIND=dp))
                 alphr = alphr + xnorm*(xnorm/real(alpha,KIND=dp))
                 tau = cmplx(alphr/beta,-alphi/beta,KIND=dp)
                 alpha = cmplx(-alphr,alphi,KIND=dp)
              end if
              alpha = la_zladiv(cmplx(one,KIND=dp),alpha)
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero (or two or whatever makes a nonnegative real number for beta).
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 alphr = real(savealpha,KIND=dp)
                 alphi = aimag(savealpha)
                 if (alphi == zero) then
                    if (alphr >= zero) then
                       tau = zero
                    else
                       tau = two
                       do j = 1,n - 1
                          x(1 + (j - 1)*incx) = zero
                       end do
                       beta = real(-savealpha,KIND=dp)
                    end if
                 else
                    xnorm = la_dlapy2(alphr,alphi)
                    tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=dp)
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    beta = xnorm
                 end if
              else
                 ! this is the general case.
                 call la_zscal(n - 1,alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_zlarfgp
#ifdef LA_WITH_XDP
     !> YLARFGP: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is real and non-negative, and
     !> x is an (n-1)-element complex vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.

     subroutine la_ylarfgp(n,alpha,x,incx,tau)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(xdp),intent(inout) :: alpha
           complex(xdp),intent(out) :: tau
           ! Array Arguments
           complex(xdp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(xdp) :: alphi,alphr,beta,bignum,smlnum,xnorm
           complex(xdp) :: savealpha
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_xynrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=xdp)
           alphi = aimag(alpha)
           if (xnorm == zero) then
              ! h  =  [1-alpha/abs(alpha) 0; 0 i], sign chosen so alpha >= 0.
              if (alphi == zero) then
                 if (alphr >= zero) then
                    ! when tau.eq.zero, the vector is special-cased to be
                    ! all zeros in the application routines.  we do not need
                    ! to clear it.
                    tau = zero
                 else
                    ! however, the application routines rely on explicit
                    ! zero checks when tau.ne.zero, and we must clear x.
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    alpha = -alpha
                 end if
              else
                 ! only "reflecting" the diagonal entry to be real and non-negative.
                 xnorm = la_xlapy2(alphr,alphi)
                 tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=xdp)
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = zero
                 end do
                 alpha = xnorm
              end if
           else
              ! general case
              beta = sign(la_xlapy3(alphr,alphi,xnorm),alphr)
              smlnum = la_xlamch('S')/la_xlamch('E')
              bignum = one/smlnum
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_yxscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alphi = alphi*bignum
                 alphr = alphr*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_xynrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=xdp)
                 beta = sign(la_xlapy3(alphr,alphi,xnorm),alphr)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alphr = alphi*(alphi/real(alpha,KIND=xdp))
                 alphr = alphr + xnorm*(xnorm/real(alpha,KIND=xdp))
                 tau = cmplx(alphr/beta,-alphi/beta,KIND=xdp)
                 alpha = cmplx(-alphr,alphi,KIND=xdp)
              end if
              alpha = la_yladiv(cmplx(one,KIND=xdp),alpha)
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero (or two or whatever makes a nonnegative real number for beta).
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 alphr = real(savealpha,KIND=xdp)
                 alphi = aimag(savealpha)
                 if (alphi == zero) then
                    if (alphr >= zero) then
                       tau = zero
                    else
                       tau = two
                       do j = 1,n - 1
                          x(1 + (j - 1)*incx) = zero
                       end do
                       beta = real(-savealpha,KIND=xdp)
                    end if
                 else
                    xnorm = la_xlapy2(alphr,alphi)
                    tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=xdp)
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    beta = xnorm
                 end if
              else
                 ! this is the general case.
                 call la_yscal(n - 1,alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_ylarfgp
#endif
#ifdef LA_WITH_QP
     !> WLARFGP: generates a complex elementary reflector H of order n, such
     !> that
     !> H**H * ( alpha ) = ( beta ),   H**H * H = I.
     !> (   x   )   (   0  )
     !> where alpha and beta are scalars, beta is real and non-negative, and
     !> x is an (n-1)-element complex vector.  H is represented in the form
     !> H = I - tau * ( 1 ) * ( 1 v**H ) ,
     !> ( v )
     !> where tau is a complex scalar and v is a complex (n-1)-element
     !> vector. Note that H is not hermitian.
     !> If the elements of x are all zero and alpha is real, then tau = 0
     !> and H is taken to be the unit matrix.

     subroutine la_wlarfgp(n,alpha,x,incx,tau)
        use la_constants_qp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           complex(qp),intent(inout) :: alpha
           complex(qp),intent(out) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,knt
           real(qp) :: alphi,alphr,beta,bignum,smlnum,xnorm
           complex(qp) :: savealpha
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag,sign
           ! Executable Statements
           if (n <= 0) then
              tau = zero
              return
           end if
           xnorm = la_qwnrm2(n - 1,x,incx)
           alphr = real(alpha,KIND=qp)
           alphi = aimag(alpha)
           if (xnorm == zero) then
              ! h  =  [1-alpha/abs(alpha) 0; 0 i], sign chosen so alpha >= 0.
              if (alphi == zero) then
                 if (alphr >= zero) then
                    ! when tau.eq.zero, the vector is special-cased to be
                    ! all zeros in the application routines.  we do not need
                    ! to clear it.
                    tau = zero
                 else
                    ! however, the application routines rely on explicit
                    ! zero checks when tau.ne.zero, and we must clear x.
                    tau = two
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    alpha = -alpha
                 end if
              else
                 ! only "reflecting" the diagonal entry to be real and non-negative.
                 xnorm = la_qlapy2(alphr,alphi)
                 tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=qp)
                 do j = 1,n - 1
                    x(1 + (j - 1)*incx) = zero
                 end do
                 alpha = xnorm
              end if
           else
              ! general case
              beta = sign(la_qlapy3(alphr,alphi,xnorm),alphr)
              smlnum = la_qlamch('S')/la_qlamch('E')
              bignum = one/smlnum
              knt = 0
              if (abs(beta) < smlnum) then
                 ! xnorm, beta may be inaccurate; scale x and recompute them
                 10 continue
                 knt = knt + 1
                 call la_wqscal(n - 1,bignum,x,incx)
                 beta = beta*bignum
                 alphi = alphi*bignum
                 alphr = alphr*bignum
                 if ((abs(beta) < smlnum) .and. (knt < 20)) go to 10
                 ! new beta is at most 1, at least smlnum
                 xnorm = la_qwnrm2(n - 1,x,incx)
                 alpha = cmplx(alphr,alphi,KIND=qp)
                 beta = sign(la_qlapy3(alphr,alphi,xnorm),alphr)
              end if
              savealpha = alpha
              alpha = alpha + beta
              if (beta < zero) then
                 beta = -beta
                 tau = -alpha/beta
              else
                 alphr = alphi*(alphi/real(alpha,KIND=qp))
                 alphr = alphr + xnorm*(xnorm/real(alpha,KIND=qp))
                 tau = cmplx(alphr/beta,-alphi/beta,KIND=qp)
                 alpha = cmplx(-alphr,alphi,KIND=qp)
              end if
              alpha = la_wladiv(cmplx(one,KIND=qp),alpha)
              if (abs(tau) <= smlnum) then
                 ! in the case where the computed tau ends up being a denormalized number,
                 ! it loses relative accuracy. this is a big problem. solution: flush tau
                 ! to zero (or two or whatever makes a nonnegative real number for beta).
                 ! (bug report provided by pat quillen from mathworks on jul 29, 2009.)
                 ! (thanks pat. thanks mathworks.)
                 alphr = real(savealpha,KIND=qp)
                 alphi = aimag(savealpha)
                 if (alphi == zero) then
                    if (alphr >= zero) then
                       tau = zero
                    else
                       tau = two
                       do j = 1,n - 1
                          x(1 + (j - 1)*incx) = zero
                       end do
                       beta = real(-savealpha,KIND=qp)
                    end if
                 else
                    xnorm = la_qlapy2(alphr,alphi)
                    tau = cmplx(one - alphr/xnorm,-alphi/xnorm,KIND=qp)
                    do j = 1,n - 1
                       x(1 + (j - 1)*incx) = zero
                    end do
                    beta = xnorm
                 end if
              else
                 ! this is the general case.
                 call la_wscal(n - 1,alpha,x,incx)
              end if
              ! if beta is subnormal, it may lose relative accuracy
              do j = 1,knt
                 beta = beta*smlnum
              end do
              alpha = beta
           end if
           return
     end subroutine la_wlarfgp
#endif

     !> CLARFT: forms the triangular factor T of a complex block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V

     pure subroutine la_clarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(sp),intent(out) :: t(ldt,*)
           complex(sp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(prevlastv,i)
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*conjg(v(i,j))
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**h * v(i:j,i)
                       call la_cgemv('CONJUGATE TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1), &
                                 ldv,v(i + 1,i),1,cone,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**h
                       call la_cgemm('N','C',i - 1,1,j - i,-tau(i),v(1,i + 1),ldv,v(i, &
                                  i + 1),ldv,cone,t(1,i),ldt)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_ctrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                  end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*conjg(v(n - k + i,j))
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**h * v(j:n-k+i,i)
                          call la_cgemv('CONJUGATE TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j, &
                                    i + 1),ldv,v(j,i),1,cone,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**h
                          call la_cgemm('N','C',k - i,1,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,cone,t(i + 1,i),ldt)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_ctrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_clarft
     !> ZLARFT: forms the triangular factor T of a complex block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V

     pure subroutine la_zlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(dp),intent(out) :: t(ldt,*)
           complex(dp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(prevlastv,i)
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*conjg(v(i,j))
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**h * v(i:j,i)
                       call la_zgemv('CONJUGATE TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1), &
                                 ldv,v(i + 1,i),1,cone,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**h
                       call la_zgemm('N','C',i - 1,1,j - i,-tau(i),v(1,i + 1),ldv,v(i, &
                                  i + 1),ldv,cone,t(1,i),ldt)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_ztrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                  end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*conjg(v(n - k + i,j))
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**h * v(j:n-k+i,i)
                          call la_zgemv('CONJUGATE TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j, &
                                    i + 1),ldv,v(j,i),1,cone,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**h
                          call la_zgemm('N','C',k - i,1,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,cone,t(i + 1,i),ldt)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_ztrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_zlarft
#ifdef LA_WITH_XDP
     !> YLARFT: forms the triangular factor T of a complex block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V

     pure subroutine la_ylarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(xdp),intent(out) :: t(ldt,*)
           complex(xdp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(prevlastv,i)
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*conjg(v(i,j))
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**h * v(i:j,i)
                       call la_ygemv('CONJUGATE TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1), &
                                 ldv,v(i + 1,i),1,cone,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**h
                       call la_ygemm('N','C',i - 1,1,j - i,-tau(i),v(1,i + 1),ldv,v(i, &
                                  i + 1),ldv,cone,t(1,i),ldt)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_ytrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                  end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*conjg(v(n - k + i,j))
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**h * v(j:n-k+i,i)
                          call la_ygemv('CONJUGATE TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j, &
                                    i + 1),ldv,v(j,i),1,cone,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**h
                          call la_ygemm('N','C',k - i,1,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,cone,t(i + 1,i),ldt)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_ytrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_ylarft
#endif
#ifdef LA_WITH_QP
     !> WLARFT: forms the triangular factor T of a complex block reflector H
     !> of order n, which is defined as a product of k elementary reflectors.
     !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
     !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
     !> If STOREV = 'C', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th column of the array V, and
     !> H  =  I - V * T * V**H
     !> If STOREV = 'R', the vector which defines the elementary reflector
     !> H(i) is stored in the i-th row of the array V, and
     !> H  =  I - V**H * T * V

     pure subroutine la_wlarft(direct,storev,n,k,v,ldv,tau,t,ldt)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,storev
           integer(ilp),intent(in) :: k,ldt,ldv,n
           ! Array Arguments
           complex(qp),intent(out) :: t(ldt,*)
           complex(qp),intent(in) :: tau(*),v(ldv,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,prevlastv,lastv
           ! Executable Statements
           ! quick return if possible
           if (n == 0) return
           if (la_lsame(direct,'F')) then
              prevlastv = n
              do i = 1,k
                 prevlastv = max(prevlastv,i)
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = 1,i
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (la_lsame(storev,'C')) then
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(lastv,i) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*conjg(v(i,j))
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(i:j,1:i-1)**h * v(i:j,i)
                       call la_wgemv('CONJUGATE TRANSPOSE',j - i,i - 1,-tau(i),v(i + 1,1), &
                                 ldv,v(i + 1,i),1,cone,t(1,i),1)
                    else
                       ! skip any trailing zeros.
                       do lastv = n,i + 1,-1
                          if (v(i,lastv) /= czero) exit
                       end do
                       do j = 1,i - 1
                          t(j,i) = -tau(i)*v(j,i)
                       end do
                       j = min(lastv,prevlastv)
                       ! t(1:i-1,i) := - tau(i) * v(1:i-1,i:j) * v(i,i:j)**h
                       call la_wgemm('N','C',i - 1,1,j - i,-tau(i),v(1,i + 1),ldv,v(i, &
                                  i + 1),ldv,cone,t(1,i),ldt)
                    end if
                    ! t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
                    call la_wtrmv('UPPER','NO TRANSPOSE','NON-UNIT',i - 1,t,ldt,t(1,i), &
                               1)
                    t(i,i) = tau(i)
                    if (i > 1) then
                       prevlastv = max(prevlastv,lastv)
                    else
                       prevlastv = lastv
                    end if
                  end if
              end do
           else
              prevlastv = 1
              do i = k,1,-1
                 if (tau(i) == czero) then
                    ! h(i)  =  i
                    do j = i,k
                       t(j,i) = czero
                    end do
                 else
                    ! general case
                    if (i < k) then
                       if (la_lsame(storev,'C')) then
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(lastv,i) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*conjg(v(n - k + i,j))
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(j:n-k+i,i+1:k)**h * v(j:n-k+i,i)
                          call la_wgemv('CONJUGATE TRANSPOSE',n - k + i - j,k - i,-tau(i),v(j, &
                                    i + 1),ldv,v(j,i),1,cone,t(i + 1,i),1)
                       else
                          ! skip any leading zeros.
                          do lastv = 1,i - 1
                             if (v(i,lastv) /= czero) exit
                          end do
                          do j = i + 1,k
                             t(j,i) = -tau(i)*v(j,n - k + i)
                          end do
                          j = max(lastv,prevlastv)
                          ! t(i+1:k,i) = -tau(i) * v(i+1:k,j:n-k+i) * v(i,j:n-k+i)**h
                          call la_wgemm('N','C',k - i,1,n - k + i - j,-tau(i),v(i + 1,j), &
                                    ldv,v(i,j),ldv,cone,t(i + 1,i),ldt)
                       end if
                       ! t(i+1:k,i) := t(i+1:k,i+1:k) * t(i+1:k,i)
                       call la_wtrmv('LOWER','NO TRANSPOSE','NON-UNIT',k - i,t(i + 1,i + 1), &
                                 ldt,t(i + 1,i),1)
                       if (i > 1) then
                          prevlastv = min(prevlastv,lastv)
                       else
                          prevlastv = lastv
                       end if
                    end if
                    t(i,i) = tau(i)
                 end if
              end do
           end if
           return
     end subroutine la_wlarft
#endif

     !> CLARFX: applies a complex elementary reflector H to a complex m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_clarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           complex(sp),intent(in) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(in) :: v(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           complex(sp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5, &
                     v6,v7,v8,v9
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (tau == czero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_clarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              v10 = conjg(v(10))
              t10 = tau*conjg(v10)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_clarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              v10 = v(10)
              t10 = tau*conjg(v10)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
410 return
     end subroutine la_clarfx
     !> ZLARFX: applies a complex elementary reflector H to a complex m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_zlarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           complex(dp),intent(in) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(in) :: v(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           complex(dp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5, &
                     v6,v7,v8,v9
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (tau == czero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_zlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              v10 = conjg(v(10))
              t10 = tau*conjg(v10)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_zlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              v10 = v(10)
              t10 = tau*conjg(v10)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_zlarfx
#ifdef LA_WITH_XDP
     !> YLARFX: applies a complex elementary reflector H to a complex m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_ylarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           complex(xdp),intent(in) :: tau
           ! Array Arguments
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(in) :: v(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           complex(xdp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5, &
                     v6,v7,v8,v9
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (tau == czero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_ylarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              v10 = conjg(v(10))
              t10 = tau*conjg(v10)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_ylarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              v10 = v(10)
              t10 = tau*conjg(v10)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_ylarfx
#endif
#ifdef LA_WITH_QP
     !> WLARFX: applies a complex elementary reflector H to a complex m by n
     !> matrix C, from either the left or the right. H is represented in the
     !> form
     !> H = I - tau * v * v**H
     !> where tau is a complex scalar and v is a complex vector.
     !> If tau = 0, then H is taken to be the unit matrix
     !> This version uses inline code if H has order < 11.

     pure subroutine la_wlarfx(side,m,n,v,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side
           integer(ilp),intent(in) :: ldc,m,n
           complex(qp),intent(in) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(in) :: v(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j
           complex(qp) :: sum,t1,t10,t2,t3,t4,t5,t6,t7,t8,t9,v1,v10,v2,v3,v4,v5, &
                     v6,v7,v8,v9
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (tau == czero) return
           if (la_lsame(side,'L')) then
              ! form  h * c, where h has order m.
              go to(10,30,50,70,90,110,130,150,170,190) m
              ! code for general m
              call la_wlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              10 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,n
                 c(1,j) = t1*c(1,j)
              end do
              go to 410
              30 continue
              ! special code for 2 x 2 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
              end do
              go to 410
              50 continue
              ! special code for 3 x 3 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
              end do
              go to 410
              70 continue
              ! special code for 4 x 4 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
              end do
              go to 410
              90 continue
              ! special code for 5 x 5 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j)

                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
              end do
              go to 410
              110 continue
              ! special code for 6 x 6 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
              end do
              go to 410
              130 continue
              ! special code for 7 x 7 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
              end do
              go to 410
              150 continue
              ! special code for 8 x 8 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
              end do
              go to 410
              170 continue
              ! special code for 9 x 9 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
              end do
              go to 410
              190 continue
              ! special code for 10 x 10 householder
              v1 = conjg(v(1))
              t1 = tau*conjg(v1)
              v2 = conjg(v(2))
              t2 = tau*conjg(v2)
              v3 = conjg(v(3))
              t3 = tau*conjg(v3)
              v4 = conjg(v(4))
              t4 = tau*conjg(v4)
              v5 = conjg(v(5))
              t5 = tau*conjg(v5)
              v6 = conjg(v(6))
              t6 = tau*conjg(v6)
              v7 = conjg(v(7))
              t7 = tau*conjg(v7)
              v8 = conjg(v(8))
              t8 = tau*conjg(v8)
              v9 = conjg(v(9))
              t9 = tau*conjg(v9)
              v10 = conjg(v(10))
              t10 = tau*conjg(v10)
              do j = 1,n
                 sum = v1*c(1,j) + v2*c(2,j) + v3*c(3,j) + v4*c(4,j) + v5*c(5,j) + &
                           v6*c(6,j) + v7*c(7,j) + v8*c(8,j) + v9*c(9,j) + v10*c(10,j)
                 c(1,j) = c(1,j) - sum*t1
                 c(2,j) = c(2,j) - sum*t2
                 c(3,j) = c(3,j) - sum*t3
                 c(4,j) = c(4,j) - sum*t4
                 c(5,j) = c(5,j) - sum*t5
                 c(6,j) = c(6,j) - sum*t6
                 c(7,j) = c(7,j) - sum*t7
                 c(8,j) = c(8,j) - sum*t8
                 c(9,j) = c(9,j) - sum*t9
                 c(10,j) = c(10,j) - sum*t10
              end do
              go to 410
           else
              ! form  c * h, where h has order n.
              go to(210,230,250,270,290,310,330,350,370,390) n
              ! code for general n
              call la_wlarf(side,m,n,v,1,tau,c,ldc,work)
              go to 410
              210 continue
              ! special code for 1 x 1 householder
              t1 = cone - tau*v(1)*conjg(v(1))
              do j = 1,m
                 c(j,1) = t1*c(j,1)
              end do
              go to 410
              230 continue
              ! special code for 2 x 2 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
              end do
              go to 410
              250 continue
              ! special code for 3 x 3 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
              end do
              go to 410
              270 continue
              ! special code for 4 x 4 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
              end do
              go to 410
              290 continue
              ! special code for 5 x 5 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5)

                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
              end do
              go to 410
              310 continue
              ! special code for 6 x 6 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
              end do
              go to 410
              330 continue
              ! special code for 7 x 7 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
              end do
              go to 410
              350 continue
              ! special code for 8 x 8 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
              end do
              go to 410
              370 continue
              ! special code for 9 x 9 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
              end do
              go to 410
              390 continue
              ! special code for 10 x 10 householder
              v1 = v(1)
              t1 = tau*conjg(v1)
              v2 = v(2)
              t2 = tau*conjg(v2)
              v3 = v(3)
              t3 = tau*conjg(v3)
              v4 = v(4)
              t4 = tau*conjg(v4)
              v5 = v(5)
              t5 = tau*conjg(v5)
              v6 = v(6)
              t6 = tau*conjg(v6)
              v7 = v(7)
              t7 = tau*conjg(v7)
              v8 = v(8)
              t8 = tau*conjg(v8)
              v9 = v(9)
              t9 = tau*conjg(v9)
              v10 = v(10)
              t10 = tau*conjg(v10)
              do j = 1,m
                 sum = v1*c(j,1) + v2*c(j,2) + v3*c(j,3) + v4*c(j,4) + v5*c(j,5) + &
                           v6*c(j,6) + v7*c(j,7) + v8*c(j,8) + v9*c(j,9) + v10*c(j,10)
                 c(j,1) = c(j,1) - sum*t1
                 c(j,2) = c(j,2) - sum*t2
                 c(j,3) = c(j,3) - sum*t3
                 c(j,4) = c(j,4) - sum*t4
                 c(j,5) = c(j,5) - sum*t5
                 c(j,6) = c(j,6) - sum*t6
                 c(j,7) = c(j,7) - sum*t7
                 c(j,8) = c(j,8) - sum*t8
                 c(j,9) = c(j,9) - sum*t9
                 c(j,10) = c(j,10) - sum*t10
              end do
              go to 410
           end if
           410 continue
           return
     end subroutine la_wlarfx
#endif

     !> CLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n Hermitian matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_clarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_sp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           complex(sp),intent(in) :: tau
           ! Array Arguments
           complex(sp),intent(inout) :: c(ldc,*)
           complex(sp),intent(in) :: v(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           complex(sp) :: alpha
           ! Executable Statements
           if (tau == czero) return
           ! form  w:= c * v
           call la_chemv(uplo,n,cone,c,ldc,v,incv,czero,work,1)
           alpha = -chalf*tau*la_cdotc(n,work,1,v,incv)
           call la_caxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_cher2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_clarfy
     !> ZLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n Hermitian matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_zlarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_dp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           complex(dp),intent(in) :: tau
           ! Array Arguments
           complex(dp),intent(inout) :: c(ldc,*)
           complex(dp),intent(in) :: v(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           complex(dp) :: alpha
           ! Executable Statements
           if (tau == czero) return
           ! form  w:= c * v
           call la_zhemv(uplo,n,cone,c,ldc,v,incv,czero,work,1)
           alpha = -chalf*tau*la_zdotc(n,work,1,v,incv)
           call la_zaxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_zher2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_zlarfy
#ifdef LA_WITH_XDP
     !> YLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n Hermitian matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_ylarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_xdp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           complex(xdp),intent(in) :: tau
           ! Array Arguments
           complex(xdp),intent(inout) :: c(ldc,*)
           complex(xdp),intent(in) :: v(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           complex(xdp) :: alpha
           ! Executable Statements
           if (tau == czero) return
           ! form  w:= c * v
           call la_yhemv(uplo,n,cone,c,ldc,v,incv,czero,work,1)
           alpha = -chalf*tau*la_ydotc(n,work,1,v,incv)
           call la_yaxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_yher2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_ylarfy
#endif
#ifdef LA_WITH_QP
     !> WLARFY: applies an elementary reflector, or Householder matrix, H,
     !> to an n x n Hermitian matrix C, from both the left and the right.
     !> H is represented in the form
     !> H = I - tau * v * v'
     !> where  tau  is a scalar and  v  is a vector.
     !> If  tau  is  zero, then  H  is taken to be the unit matrix.

     pure subroutine la_wlarfy(uplo,n,v,incv,tau,c,ldc,work)
        use la_constants_qp
        ! -- lapack test routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incv,ldc,n
           complex(qp),intent(in) :: tau
           ! Array Arguments
           complex(qp),intent(inout) :: c(ldc,*)
           complex(qp),intent(in) :: v(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           complex(qp) :: alpha
           ! Executable Statements
           if (tau == czero) return
           ! form  w:= c * v
           call la_whemv(uplo,n,cone,c,ldc,v,incv,czero,work,1)
           alpha = -chalf*tau*la_wdotc(n,work,1,v,incv)
           call la_waxpy(n,alpha,v,incv,work,1)
           ! c := c - v * w' - w * v'
           call la_wher2(uplo,n,-tau,v,incv,work,1,c,ldc)
           return
     end subroutine la_wlarfy
#endif

end module la_lapack_householder_reflectors
