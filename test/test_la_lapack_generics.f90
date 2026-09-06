! Kind-agnostic LAPACK generics of the blocked and expert drivers
module test_la_lapack_generics
    use la_constants,only:sp,dp,qp,ilp,lk
    use la_lapack,only:geqp3,ggev3,gges3,gghd3,geevx,geesx,gesvx
    use testdrive,only:error_type,check,new_unittest,unittest_type

    implicit none(type,external)
    private

    public :: test_lapack_generics

    contains

    !> Every generic the umbrella gained, called once per kind
    subroutine test_lapack_generics(tests)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: tests(:)

        allocate (tests(0))

        call add_test(tests,new_unittest("lapack_generics_s",test_s_generics))
        call add_test(tests,new_unittest("lapack_generics_d",test_d_generics))
        call add_test(tests,new_unittest("lapack_generics_q",test_q_generics))
        call add_test(tests,new_unittest("lapack_generics_c",test_c_generics))
        call add_test(tests,new_unittest("lapack_generics_z",test_z_generics))
        call add_test(tests,new_unittest("lapack_generics_w",test_w_generics))

    end subroutine test_lapack_generics

    !> Eigenvalue selector of the generalized drivers, real(sp)
    pure logical(lk) function keep_s(alphar,alphai,beta)
        real(sp),intent(in) :: alphar,alphai,beta
        keep_s = alphar > 0.0_sp .and. alphai >= 0.0_sp .and. beta > 0.0_sp
    end function keep_s

    !> Eigenvalue selector of the Schur drivers, real(sp)
    pure logical(lk) function keep_schur_s(alphar,alphai)
        real(sp),intent(in) :: alphar,alphai
        keep_schur_s = alphar > 0.0_sp .and. alphai >= 0.0_sp
    end function keep_schur_s

    !> Eigenvalue selector of the generalized drivers, real(dp)
    pure logical(lk) function keep_d(alphar,alphai,beta)
        real(dp),intent(in) :: alphar,alphai,beta
        keep_d = alphar > 0.0_dp .and. alphai >= 0.0_dp .and. beta > 0.0_dp
    end function keep_d

    !> Eigenvalue selector of the Schur drivers, real(dp)
    pure logical(lk) function keep_schur_d(alphar,alphai)
        real(dp),intent(in) :: alphar,alphai
        keep_schur_d = alphar > 0.0_dp .and. alphai >= 0.0_dp
    end function keep_schur_d

    !> Eigenvalue selector of the generalized drivers, real(qp)
    pure logical(lk) function keep_q(alphar,alphai,beta)
        real(qp),intent(in) :: alphar,alphai,beta
        keep_q = alphar > 0.0_qp .and. alphai >= 0.0_qp .and. beta > 0.0_qp
    end function keep_q

    !> Eigenvalue selector of the Schur drivers, real(qp)
    pure logical(lk) function keep_schur_q(alphar,alphai)
        real(qp),intent(in) :: alphar,alphai
        keep_schur_q = alphar > 0.0_qp .and. alphai >= 0.0_qp
    end function keep_schur_q

    !> Eigenvalue selector of the generalized drivers, complex(sp)
    pure logical(lk) function keep_c(alpha,beta)
        complex(sp),intent(in) :: alpha,beta
        keep_c = real(alpha,sp) > 0.0_sp .and. real(beta,sp) > 0.0_sp
    end function keep_c

    !> Eigenvalue selector of the Schur drivers, complex(sp)
    pure logical(lk) function keep_schur_c(alpha)
        complex(sp),intent(in) :: alpha
        keep_schur_c = real(alpha,sp) > 0.0_sp
    end function keep_schur_c

    !> Eigenvalue selector of the generalized drivers, complex(dp)
    pure logical(lk) function keep_z(alpha,beta)
        complex(dp),intent(in) :: alpha,beta
        keep_z = real(alpha,dp) > 0.0_dp .and. real(beta,dp) > 0.0_dp
    end function keep_z

    !> Eigenvalue selector of the Schur drivers, complex(dp)
    pure logical(lk) function keep_schur_z(alpha)
        complex(dp),intent(in) :: alpha
        keep_schur_z = real(alpha,dp) > 0.0_dp
    end function keep_schur_z

    !> Eigenvalue selector of the generalized drivers, complex(qp)
    pure logical(lk) function keep_w(alpha,beta)
        complex(qp),intent(in) :: alpha,beta
        keep_w = real(alpha,qp) > 0.0_qp .and. real(beta,qp) > 0.0_qp
    end function keep_w

    !> Eigenvalue selector of the Schur drivers, complex(qp)
    pure logical(lk) function keep_schur_w(alpha)
        complex(qp),intent(in) :: alpha
        keep_schur_w = real(alpha,qp) > 0.0_qp
    end function keep_schur_w

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, real(sp)
    subroutine test_s_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n),iwork(8*n)
        logical(lk) :: bwork(n)
        character :: equed
        real(sp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        real(sp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        real(sp) :: rhs(n,1),bmat(n,1),x(n,1)
        real(sp) :: tau(n),work(lwork),alphar(n),alphai(n),beta(n),wr(n),wi(n)
        real(sp) :: scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(sp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([real(sp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_sp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_sp
        end do
        rhs(:,1) = [real(sp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,info)
        call check(error,info == 0,'geqp3 s failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 s did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alphar,alphai,beta,vl,n,vr,n,work,lwork,info)
        call check(error,info == 0,'ggev3 s failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_s,n,acopy,n,bcopy,n,sdim,alphar,alphai,beta, &
                   vsl,n,vsr,n,work,lwork,bwork,info)
        call check(error,info == 0,'gges3 s failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 s failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,wr,wi,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,iwork,info)
        call check(error,info == 0,'geevx s failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_s,'N',n,acopy,n,sdim,wr,wi,vs,n,rce,rcv, &
                   work,lwork,iwork,size(iwork,kind=ilp),bwork,info)
        call check(error,info == 0,'geesx s failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,iwork,info)
        call check(error,info == 0,'gesvx s failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx s residual too large')
        if (allocated(error)) return

    end subroutine test_s_generics

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, real(dp)
    subroutine test_d_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n),iwork(8*n)
        logical(lk) :: bwork(n)
        character :: equed
        real(dp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        real(dp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        real(dp) :: rhs(n,1),bmat(n,1),x(n,1)
        real(dp) :: tau(n),work(lwork),alphar(n),alphai(n),beta(n),wr(n),wi(n)
        real(dp) :: scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(dp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([real(dp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_dp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_dp
        end do
        rhs(:,1) = [real(dp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,info)
        call check(error,info == 0,'geqp3 d failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 d did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alphar,alphai,beta,vl,n,vr,n,work,lwork,info)
        call check(error,info == 0,'ggev3 d failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_d,n,acopy,n,bcopy,n,sdim,alphar,alphai,beta, &
                   vsl,n,vsr,n,work,lwork,bwork,info)
        call check(error,info == 0,'gges3 d failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 d failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,wr,wi,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,iwork,info)
        call check(error,info == 0,'geevx d failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_d,'N',n,acopy,n,sdim,wr,wi,vs,n,rce,rcv, &
                   work,lwork,iwork,size(iwork,kind=ilp),bwork,info)
        call check(error,info == 0,'geesx d failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,iwork,info)
        call check(error,info == 0,'gesvx d failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx d residual too large')
        if (allocated(error)) return

    end subroutine test_d_generics

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, real(qp)
    subroutine test_q_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n),iwork(8*n)
        logical(lk) :: bwork(n)
        character :: equed
        real(qp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        real(qp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        real(qp) :: rhs(n,1),bmat(n,1),x(n,1)
        real(qp) :: tau(n),work(lwork),alphar(n),alphai(n),beta(n),wr(n),wi(n)
        real(qp) :: scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(qp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([real(qp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_qp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_qp
        end do
        rhs(:,1) = [real(qp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,info)
        call check(error,info == 0,'geqp3 q failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 q did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alphar,alphai,beta,vl,n,vr,n,work,lwork,info)
        call check(error,info == 0,'ggev3 q failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_q,n,acopy,n,bcopy,n,sdim,alphar,alphai,beta, &
                   vsl,n,vsr,n,work,lwork,bwork,info)
        call check(error,info == 0,'gges3 q failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 q failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,wr,wi,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,iwork,info)
        call check(error,info == 0,'geevx q failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_q,'N',n,acopy,n,sdim,wr,wi,vs,n,rce,rcv, &
                   work,lwork,iwork,size(iwork,kind=ilp),bwork,info)
        call check(error,info == 0,'geesx q failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,iwork,info)
        call check(error,info == 0,'gesvx q failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx q residual too large')
        if (allocated(error)) return

    end subroutine test_q_generics

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, complex(sp)
    subroutine test_c_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(sp),parameter :: tol = 100*sqrt(epsilon(0.0_sp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n)
        logical(lk) :: bwork(n)
        character :: equed
        complex(sp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        complex(sp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        complex(sp) :: rhs(n,1),bmat(n,1),x(n,1)
        complex(sp) :: tau(n),work(lwork),alpha(n),beta(n),w(n)
        real(sp) :: rwork(8*n),scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(sp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([complex(sp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_sp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_sp
        end do
        rhs(:,1) = [complex(sp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,rwork,info)
        call check(error,info == 0,'geqp3 c failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 c did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alpha,beta,vl,n,vr,n,work,lwork,rwork,info)
        call check(error,info == 0,'ggev3 c failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_c,n,acopy,n,bcopy,n,sdim,alpha,beta, &
                   vsl,n,vsr,n,work,lwork,rwork,bwork,info)
        call check(error,info == 0,'gges3 c failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 c failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,w,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,rwork,info)
        call check(error,info == 0,'geevx c failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_c,'N',n,acopy,n,sdim,w,vs,n,rce,rcv, &
                   work,lwork,rwork,bwork,info)
        call check(error,info == 0,'geesx c failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,rwork,info)
        call check(error,info == 0,'gesvx c failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx c residual too large')
        if (allocated(error)) return

    end subroutine test_c_generics

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, complex(dp)
    subroutine test_z_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(dp),parameter :: tol = 100*sqrt(epsilon(0.0_dp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n)
        logical(lk) :: bwork(n)
        character :: equed
        complex(dp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        complex(dp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        complex(dp) :: rhs(n,1),bmat(n,1),x(n,1)
        complex(dp) :: tau(n),work(lwork),alpha(n),beta(n),w(n)
        real(dp) :: rwork(8*n),scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(dp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([complex(dp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_dp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_dp
        end do
        rhs(:,1) = [complex(dp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,rwork,info)
        call check(error,info == 0,'geqp3 z failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 z did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alpha,beta,vl,n,vr,n,work,lwork,rwork,info)
        call check(error,info == 0,'ggev3 z failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_z,n,acopy,n,bcopy,n,sdim,alpha,beta, &
                   vsl,n,vsr,n,work,lwork,rwork,bwork,info)
        call check(error,info == 0,'gges3 z failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 z failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,w,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,rwork,info)
        call check(error,info == 0,'geevx z failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_z,'N',n,acopy,n,sdim,w,vs,n,rce,rcv, &
                   work,lwork,rwork,bwork,info)
        call check(error,info == 0,'geesx z failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,rwork,info)
        call check(error,info == 0,'gesvx z failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx z residual too large')
        if (allocated(error)) return

    end subroutine test_z_generics

    !> geqp3, ggev3, gges3, gghd3, geevx, geesx and gesvx, complex(qp)
    subroutine test_w_generics(error)
        type(error_type),allocatable,intent(out) :: error

        integer(ilp),parameter :: n = 3_ilp,lwork = 512_ilp
        real(qp),parameter :: tol = 100*sqrt(epsilon(0.0_qp))

        integer(ilp) :: i,info,sdim,ilo,ihi,jpvt(n),ipiv(n)
        logical(lk) :: bwork(n)
        character :: equed
        complex(qp) :: a(n,n),eye(n,n),acopy(n,n),bcopy(n,n),af(n,n)
        complex(qp) :: q(n,n),z(n,n),vl(n,n),vr(n,n),vs(n,n),vsl(n,n),vsr(n,n)
        complex(qp) :: rhs(n,1),bmat(n,1),x(n,1)
        complex(qp) :: tau(n),work(lwork),alpha(n),beta(n),w(n)
        real(qp) :: rwork(8*n),scal(n),rconde(n),rcondv(n),rowsc(n),colsc(n),ferr(1),berr(1)
        real(qp) :: abnrm,rcond,rce,rcv

        a = transpose(reshape([complex(qp) :: 4,1,0, &
                                        1,3,1, &
                                        0,1,2], [n,n]))
        eye = 0.0_qp
        do concurrent(i=1:n)
          eye(i,i) = 1.0_qp
        end do
        rhs(:,1) = [complex(qp) :: 1,2,3]

        !> Pivoted QR factorization
        acopy = a
        jpvt = 0_ilp
        call geqp3(n,n,acopy,n,jpvt,tau,work,lwork,rwork,info)
        call check(error,info == 0,'geqp3 w failed')
        if (allocated(error)) return
        call check(error,all(jpvt >= 1_ilp) .and. all(jpvt <= n) .and. sum(jpvt) == n*(n + 1)/2, &
                   'geqp3 w did not return a permutation')
        if (allocated(error)) return

        !> Generalized eigenvalues, blocked
        acopy = a
        bcopy = eye
        call ggev3('N','V',n,acopy,n,bcopy,n,alpha,beta,vl,n,vr,n,work,lwork,rwork,info)
        call check(error,info == 0,'ggev3 w failed')
        if (allocated(error)) return

        !> Generalized Schur form, blocked
        acopy = a
        bcopy = eye
        call gges3('N','N','N',keep_w,n,acopy,n,bcopy,n,sdim,alpha,beta, &
                   vsl,n,vsr,n,work,lwork,rwork,bwork,info)
        call check(error,info == 0,'gges3 w failed')
        if (allocated(error)) return

        !> Generalized Hessenberg form, blocked
        acopy = a
        bcopy = eye
        call gghd3('N','N',n,1_ilp,n,acopy,n,bcopy,n,q,n,z,n,work,lwork,info)
        call check(error,info == 0,'gghd3 w failed')
        if (allocated(error)) return

        !> Eigenvalues with condition numbers
        acopy = a
        call geevx('B','V','V','B',n,acopy,n,w,vl,n,vr,n,ilo,ihi,scal,abnrm, &
                   rconde,rcondv,work,lwork,rwork,info)
        call check(error,info == 0,'geevx w failed')
        if (allocated(error)) return

        !> Schur form with condition numbers
        acopy = a
        call geesx('N','N',keep_schur_w,'N',n,acopy,n,sdim,w,vs,n,rce,rcv, &
                   work,lwork,rwork,bwork,info)
        call check(error,info == 0,'geesx w failed')
        if (allocated(error)) return

        !> Expert linear solve
        acopy = a
        bmat = rhs
        equed = 'N'
        call gesvx('N','N',n,1_ilp,acopy,n,af,n,ipiv,equed,rowsc,colsc,bmat,n,x,n, &
                   rcond,ferr,berr,work,rwork,info)
        call check(error,info == 0,'gesvx w failed')
        if (allocated(error)) return
        call check(error,all(abs(matmul(a,x(:,1)) - rhs(:,1)) < tol),'gesvx w residual too large')
        if (allocated(error)) return

    end subroutine test_w_generics

    ! gcc-15 bugfix utility
    subroutine add_test(tests,new_test)
        type(unittest_type),allocatable,intent(inout) :: tests(:)
        type(unittest_type),intent(in) :: new_test

        integer :: n
        type(unittest_type),allocatable :: new_tests(:)

        if (allocated(tests)) then
            n = size(tests)
        else
            n = 0
        end if

        allocate (new_tests(n + 1))
        if (n > 0) new_tests(1:n) = tests(1:n)
                 new_tests(1 + n) = new_test
        call move_alloc(from=new_tests,to=tests)

    end subroutine add_test

end module test_la_lapack_generics

program test_lapack_generics_program
     use,intrinsic :: iso_fortran_env,only:error_unit
     use testdrive,only:run_testsuite,new_testsuite,testsuite_type
     use test_la_lapack_generics,only:test_lapack_generics
     implicit none
     integer :: stat,is
     type(testsuite_type),allocatable :: testsuites(:)
     character(len=*),parameter :: fmt = '("#", *(1x, a))'

     stat = 0

     testsuites = [ &
         new_testsuite("la_lapack_generics",test_lapack_generics) &
         ]

     do is = 1,size(testsuites)
         write (error_unit,fmt) "Testing:",testsuites(is)%name
         call run_testsuite(testsuites(is)%collect,error_unit,stat)
     end do

     if (stat > 0) then
         write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
         error stop
     end if
end program test_lapack_generics_program
