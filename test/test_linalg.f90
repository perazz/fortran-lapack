
module test_linalg
    use testdrive,only:new_unittest,unittest_type,error_type,check
    use la_constants,only:sp,dp,qp
    use linear_algebra,only:diag,eye,trace,outer_product,cross_product,kronecker_product,hermitian

    implicit none

    real(sp),parameter :: sptol = 1000*epsilon(1._sp)
    real(dp),parameter :: dptol = 1000*epsilon(1._dp)
    real(qp),parameter :: qptol = 1000*epsilon(1._qp)

contains

    !> Collect all exported unit tests
    subroutine collect_linalg(testsuite)
        !> Collection of tests
        type(unittest_type),allocatable,intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("eye",test_eye), &
            new_unittest("diag_rsp",test_diag_rsp), &
            new_unittest("diag_rsp_k",test_diag_rsp_k), &
            new_unittest("diag_rdp",test_diag_rdp), &
            new_unittest("diag_rqp",test_diag_rqp), &
            new_unittest("diag_csp",test_diag_csp), &
            new_unittest("diag_cdp",test_diag_cdp), &
            new_unittest("diag_cqp",test_diag_cqp), &
            new_unittest("trace_rsp",test_trace_rsp), &
            new_unittest("trace_rsp_nonsquare",test_trace_rsp_nonsquare), &
            new_unittest("trace_rdp",test_trace_rdp), &
            new_unittest("trace_rdp_nonsquare",test_trace_rdp_nonsquare), &
            new_unittest("trace_rqp",test_trace_rqp), &
            new_unittest("trace_csp",test_trace_csp), &
            new_unittest("trace_cdp",test_trace_cdp), &
            new_unittest("trace_cqp",test_trace_cqp), &
            new_unittest("kronecker_product_rsp",test_kronecker_product_rsp), &
            new_unittest("kronecker_product_rdp",test_kronecker_product_rdp), &
            new_unittest("kronecker_product_rqp",test_kronecker_product_rqp), &
            new_unittest("kronecker_product_csp",test_kronecker_product_csp), &
            new_unittest("kronecker_product_cdp",test_kronecker_product_cdp), &
            new_unittest("kronecker_product_cqp",test_kronecker_product_cqp), &
            new_unittest("hermitian_csp",test_hermitian_csp), &
            new_unittest("hermitian_cdp",test_hermitian_cdp), &
            new_unittest("hermitian_cqp",test_hermitian_cqp), &
            new_unittest("outer_product_rsp",test_outer_product_rsp), &
            new_unittest("outer_product_rdp",test_outer_product_rdp), &
            new_unittest("outer_product_rqp",test_outer_product_rqp), &
            new_unittest("outer_product_csp",test_outer_product_csp), &
            new_unittest("outer_product_cdp",test_outer_product_cdp), &
            new_unittest("outer_product_cqp",test_outer_product_cqp), &
            new_unittest("cross_product_rsp",test_cross_product_rsp), &
            new_unittest("cross_product_rdp",test_cross_product_rdp), &
            new_unittest("cross_product_rqp",test_cross_product_rqp), &
            new_unittest("cross_product_csp",test_cross_product_csp), &
            new_unittest("cross_product_cdp",test_cross_product_cdp), &
            new_unittest("cross_product_cqp",test_cross_product_cqp) &
            ]

    end subroutine collect_linalg

    subroutine test_eye(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        real(sp),allocatable :: rye(:,:)
        complex(sp) :: cye(7,7)
        integer :: i

        call check(error,all(eye(3,3) == diag([(1.0_dp,i=1,3)])), &
            "all(eye(3,3) == diag([(1.0_dp,i=1,3)])) failed.")
        if (allocated(error)) return

        rye = real(eye(3,4),sp)
        call check(error,sum(abs(rye(:,1:3) - diag([(1.0_sp,i=1,3)]))) < sptol, &
            "sum(abs(rye(:,1:3) - diag([(1.0_sp,i=1,3)]))) < sptol failed")
        if (allocated(error)) return

        call check(error,all(eye(5) == diag([(1.0_dp,i=1,5)])), &
            "all(eye(5) == diag([(1.0_dp,i=1,5)] failed.")
        if (allocated(error)) return

        rye = real(eye(6),sp)
        call check(error,sum(rye - diag([(1.0_sp,i=1,6)])) < sptol, &
            "sum(rye - diag([(1.0_sp,i=1,6)])) < sptol failed.")
        if (allocated(error)) return

        cye = real(eye(7),sp)
        call check(error,abs(trace(cye) - cmplx(7.0_sp,0.0_sp,kind=sp)) < sptol, &
            "abs(trace(cye) - cmplx(7.0_sp,0.0_sp,kind=sp)) < sptol failed.")

    end subroutine test_eye

    subroutine test_diag_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")
        if (allocated(error)) return

        call check(error,all(diag(3*a) == 3*v), &
            "all(diag(3*a) == 3*v) failed.")
    end subroutine test_diag_rsp

    subroutine test_diag_rsp_k(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 4
        real(sp) :: a(n,n),b(n,n)
        integer :: i,j

        a = diag([(1._sp,i=1,n - 1)],-1)

        b = reshape([((merge(1,0,i == j + 1),i=1,n),j=1,n)], [n,n])

        call check(error,all(a == b), &
            "all(a == b) failed.")
        if (allocated(error)) return

        call check(error,sum(diag(a,-1)) - (n - 1) < sptol, &
            "sum(diag(a,-1)) - (n-1) < sptol failed.")
        if (allocated(error)) return

        call check(error,all(a == transpose(diag([(1._sp,i=1,n - 1)],1))), &
            "all(a == transpose(diag([(1._sp,i=1,n-1)],1))) failed")
        if (allocated(error)) return

        call random_number(a)
        do i = 1,n
            call check(error,size(diag(a,i)) == n - i, &
                "size(diag(a,i)) == n-i failed.")
            if (allocated(error)) return
        end do
        call check(error,size(diag(a,n + 1)) == 0, &
            "size(diag(a,n+1)) == 0 failed.")
    end subroutine test_diag_rsp_k

    subroutine test_diag_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")
        if (allocated(error)) return

        call check(error,all(diag(3*a) == 3*v), &
            "all(diag(3*a) == 3*v) failed.")
    end subroutine test_diag_rdp

    subroutine test_diag_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp) :: v(n),a(n,n),b(n,n)
        integer :: i,j

        v = [(i,i=1,n)]
        a = diag(v)
        b = reshape([((merge(i,0,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")
        if (allocated(error)) return

        call check(error,all(diag(3*a) == 3*v), &
            "all(diag(3*a) == 3*v) failed.")
    end subroutine test_diag_rqp

    subroutine test_diag_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(sp) :: a(n,n),b(n,n)
        complex(sp),parameter :: i_ = cmplx(0,1,kind=sp)
        integer :: i,j

        a = diag([(cmplx(i,0,kind=sp),i=1,n)]) + diag([(i_,i=1,n)])
        b = reshape([((merge(i + 1*i_,0*i_,i == j),i=1,n),j=1,n)], [n,n])
        call check(error,all(a == b), &
            "all(a == b) failed.")
        if (allocated(error)) return

        call check(error,all(abs(real(diag(a)) - [(i,i=1,n)]) < sptol), &
            "all(abs(real(diag(a)) - [(i,i=1,n)]) < sptol)")
        if (allocated(error)) return
        call check(error,all(abs(aimag(diag(a)) - [(1,i=1,n)]) < sptol), &
            "all(abs(aimag(diag(a)) - [(1,i=1,n)]) < sptol)")
    end subroutine test_diag_csp

    subroutine test_diag_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(dp) :: a(n,n)
        complex(dp),parameter :: i_ = cmplx(0,1,kind=dp)

        a = diag([i_],-2) + diag([i_],2)
        call check(error,a(3,1) == i_ .and. a(1,3) == i_, &
            "a(3,1) == i_ .and. a(1,3) == i_ failed.")
    end subroutine test_diag_cdp

    subroutine test_diag_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(qp) :: a(n,n)
        complex(qp),parameter :: i_ = cmplx(0,1,kind=qp)

        a = diag([i_,i_],-1) + diag([i_,i_],1)
        call check(error,all(diag(a,-1) == i_) .and. all(diag(a,1) == i_), &
            "all(diag(a,-1) == i_) .and. all(diag(a,1) == i_) failed.")
    end subroutine test_diag_cqp

    subroutine test_trace_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 5
        real(sp) :: a(n,n)
        integer :: i

        a = reshape([(i,i=1,n**2)], [n,n])
        call check(error,abs(trace(a) - sum(diag(a))) < sptol, &
            "abs(trace(a) - sum(diag(a))) < sptol failed.")
    end subroutine test_trace_rsp

    subroutine test_trace_rsp_nonsquare(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 4
        real(sp) :: a(n,n + 1),ans
        integer :: i

        ! 1   5   9  13  17
        ! 2   6  10  14  18
        ! 3   7  11  15  19
        ! 4   8  12  16  20
        a = reshape([(i,i=1,n*(n + 1))], [n,n + 1])
        ans = sum([1._sp,6._sp,11._sp,16._sp])

        call check(error,abs(trace(a) - ans) < sptol, &
            "abs(trace(a) - ans) < sptol failed.")
    end subroutine test_trace_rsp_nonsquare

    subroutine test_trace_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 4
        real(dp) :: a(n,n)
        integer :: i

        a = reshape([(i,i=1,n**2)], [n,n])
        call check(error,abs(trace(a) - sum(diag(a))) < dptol, &
            "abs(trace(a) - sum(diag(a))) < dptol failed.")
    end subroutine test_trace_rdp

    subroutine test_trace_rdp_nonsquare(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 4
        real(dp) :: a(n,n - 1),ans
        integer :: i

        !  1   25   81
        !  4   36  100
        !  9   49  121
        ! 16   64  144
        a = reshape([(i**2,i=1,n*(n - 1))], [n,n - 1])
        ans = sum([1._dp,36._dp,121._dp])

        call check(error,abs(trace(a) - ans) < dptol, &
            "abs(trace(a) - ans) < dptol failed.")
    end subroutine test_trace_rdp_nonsquare

    subroutine test_trace_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp) :: a(n,n)
        integer :: i

        a = reshape([(i,i=1,n**2)], [n,n])
        call check(error,abs(trace(a) - sum(diag(a))) < qptol, &
            "abs(trace(a) - sum(diag(a))) < qptol failed.")
    end subroutine test_trace_rqp

    subroutine test_trace_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 5
        real(sp) :: re(n,n),im(n,n)
        complex(sp) :: a(n,n),b(n,n)
        complex(sp),parameter :: i_ = cmplx(0,1,kind=sp)

        call random_number(re)
        call random_number(im)
        a = re + im*i_

        call random_number(re)
        call random_number(im)
        b = re + im*i_

        ! tr(A + B) = tr(A) + tr(B)
        call check(error,abs(trace(a + b) - (trace(a) + trace(b))) < sptol, &
            "abs(trace(a+b) - (trace(a) + trace(b))) < sptol failed.")
    end subroutine test_trace_csp

    subroutine test_trace_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(dp) :: a(n,n),ans
        complex(dp),parameter :: i_ = cmplx(0,1,kind=dp)
        integer :: j

        a = reshape([(j + (n**2 - (j - 1))*i_,j=1,n**2)], [n,n])
        ans = cmplx(15,15,kind=dp) !(1 + 5 + 9) + (9 + 5 + 1)i

        call check(error,abs(trace(a) - ans) < dptol, &
            "abs(trace(a) - ans) < dptol failed.")
    end subroutine test_trace_cdp

    subroutine test_trace_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(qp) :: a(n,n)
        complex(qp),parameter :: i_ = cmplx(0,1,kind=qp)

        a = 3*eye(n) + 4*eye(n)*i_ ! pythagorean triple
        call check(error,abs(trace(a)) - 3*5.0_qp < qptol, &
            "abs(trace(a)) - 3*5.0_qp < qptol failed.")
    end subroutine test_trace_cqp

    subroutine test_kronecker_product_rsp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      real(sp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      real(sp),parameter :: tol = 1.e-6

      real(sp) :: A(m1,n1),B(m2,n2)
      real(sp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_rsp
    subroutine test_kronecker_product_rdp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      real(dp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      real(dp),parameter :: tol = 1.e-6

      real(dp) :: A(m1,n1),B(m2,n2)
      real(dp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_rdp
    subroutine test_kronecker_product_rqp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      real(qp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      real(qp),parameter :: tol = 1.e-6

      real(qp) :: A(m1,n1),B(m2,n2)
      real(qp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_rqp
    subroutine test_kronecker_product_csp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      complex(sp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      complex(sp),parameter :: tol = 1.e-6

      complex(sp) :: A(m1,n1),B(m2,n2)
      complex(sp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_csp
    subroutine test_kronecker_product_cdp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      complex(dp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      complex(dp),parameter :: tol = 1.e-6

      complex(dp) :: A(m1,n1),B(m2,n2)
      complex(dp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_cdp
    subroutine test_kronecker_product_cqp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m1 = 1,n1 = 2,m2 = 2,n2 = 3
      complex(qp),dimension(m1*m2,n1*n2),parameter :: expected &
           = transpose(reshape([1,2,3,2,4,6,2,4,6,4,8,12], [m2*n2,m1*n1]))
      complex(qp),parameter :: tol = 1.e-6

      complex(qp) :: A(m1,n1),B(m2,n2)
      complex(qp) :: C(m1*m2,n1*n2),diff(m1*m2,n1*n2)

      integer :: i,j

      do j = 1,n1
         do i = 1,m1
            A(i,j) = i*j ! A = [1, 2]
         end do
      end do

      do j = 1,n2
         do i = 1,m2
            B(i,j) = i*j ! B = [[1, 2, 3], [2, 4, 6]]
         end do
      end do

      C = kronecker_product(A,B)

      diff = C - expected

      call check(error,all(abs(diff) <= abs(tol)),"all(abs(diff) .le. abs(tol)) failed")
      ! Expected: C = [1*B, 2*B] = [[1,2,3, 2,4,6], [2,4,6, 4, 8, 12]]

    end subroutine test_kronecker_product_cqp

    subroutine test_hermitian_csp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m = 2,n = 3
      complex(sp),dimension(m,n) :: A
      complex(sp),dimension(n,m) :: AT,expected,diff
      real(sp),parameter :: tol = 1.e-6_sp

      integer :: i,j
      
      do concurrent(i=1:m,j=1:n)
        A(i,j) = cmplx(i,-j,kind=sp)
        expected(j,i) = cmplx(i,+j,kind=sp)
      end do

      AT = hermitian(A)

      diff = AT - expected

      call check(error,all(abs(diff) < abs(tol)),"hermitian: all(abs(diff) < abs(tol)) failed")

    end subroutine test_hermitian_csp
    subroutine test_hermitian_cdp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m = 2,n = 3
      complex(dp),dimension(m,n) :: A
      complex(dp),dimension(n,m) :: AT,expected,diff
      real(dp),parameter :: tol = 1.e-6_dp

      integer :: i,j
      
      do concurrent(i=1:m,j=1:n)
        A(i,j) = cmplx(i,-j,kind=dp)
        expected(j,i) = cmplx(i,+j,kind=dp)
      end do

      AT = hermitian(A)

      diff = AT - expected

      call check(error,all(abs(diff) < abs(tol)),"hermitian: all(abs(diff) < abs(tol)) failed")

    end subroutine test_hermitian_cdp
    subroutine test_hermitian_cqp(error)
      !> Error handling
      type(error_type),allocatable,intent(out) :: error
      integer,parameter :: m = 2,n = 3
      complex(qp),dimension(m,n) :: A
      complex(qp),dimension(n,m) :: AT,expected,diff
      real(qp),parameter :: tol = 1.e-6_qp

      integer :: i,j
      
      do concurrent(i=1:m,j=1:n)
        A(i,j) = cmplx(i,-j,kind=qp)
        expected(j,i) = cmplx(i,+j,kind=qp)
      end do

      AT = hermitian(A)

      diff = AT - expected

      call check(error,all(abs(diff) < abs(tol)),"hermitian: all(abs(diff) < abs(tol)) failed")

    end subroutine test_hermitian_cqp

    subroutine test_outer_product_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        real(sp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [1.,2.]
        v = [1.,3.]
        expected = reshape([1.,2.,3.,6.], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < sptol), &
            "all(abs(diff) < sptol) failed.")
    end subroutine test_outer_product_rsp

    subroutine test_outer_product_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        real(dp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [1.,2.]
        v = [1.,3.]
        expected = reshape([1.,2.,3.,6.], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < dptol), &
            "all(abs(diff) < dptol) failed.")
    end subroutine test_outer_product_rdp

    subroutine test_outer_product_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        real(qp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [1.,2.]
        v = [1.,3.]
        expected = reshape([1.,2.,3.,6.], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < qptol), &
            "all(abs(diff) < qptol) failed.")
    end subroutine test_outer_product_rqp

    subroutine test_outer_product_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        complex(sp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [cmplx(1.,1.),cmplx(2.,0.)]
        v = [cmplx(1.,0.),cmplx(3.,1.)]
        expected = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(2.,4.),cmplx(6.,2.)], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < sptol), &
            "all(abs(diff) < sptol) failed.")
    end subroutine test_outer_product_csp

    subroutine test_outer_product_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        complex(dp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [cmplx(1.,1.),cmplx(2.,0.)]
        v = [cmplx(1.,0.),cmplx(3.,1.)]
        expected = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(2.,4.),cmplx(6.,2.)], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < dptol), &
            "all(abs(diff) < dptol) failed.")
    end subroutine test_outer_product_cdp

    subroutine test_outer_product_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 2
        complex(qp) :: u(n),v(n),expected(n,n),diff(n,n)

        u = [cmplx(1.,1.),cmplx(2.,0.)]
        v = [cmplx(1.,0.),cmplx(3.,1.)]
        expected = reshape([cmplx(1.,1.),cmplx(2.,0.),cmplx(2.,4.),cmplx(6.,2.)], [n,n])
        diff = expected - outer_product(u,v)
        call check(error,all(abs(diff) < qptol), &
            "all(abs(diff) < qptol) failed.")
    end subroutine test_outer_product_cqp

    subroutine test_cross_product_rsp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(sp) :: u(n),v(n),expected(n),diff(n)
        u = [1.1_sp,2.5_sp,2.4_sp]
        v = [0.5_sp,1.5_sp,2.5_sp]
        expected = [2.65_sp,-1.55_sp,0.4_sp]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < sptol), &
             "all(abs(cross_product(u,v)-expected)) < sptol failed.")
    end subroutine test_cross_product_rsp

    subroutine test_cross_product_rdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(dp) :: u(n),v(n),expected(n),diff(n)
        u = [1.1_dp,2.5_dp,2.4_dp]
        v = [0.5_dp,1.5_dp,2.5_dp]
        expected = [2.65_dp,-1.55_dp,0.4_dp]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < dptol), &
             "all(abs(cross_product(u,v)-expected)) < dptol failed.")
    end subroutine test_cross_product_rdp

    subroutine test_cross_product_rqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        real(qp) :: u(n),v(n),expected(n),diff(n)
        u = [1.1_qp,2.5_qp,2.4_qp]
        v = [0.5_qp,1.5_qp,2.5_qp]
        expected = [2.65_qp,-1.55_qp,0.4_qp]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < qptol), &
             "all(abs(cross_product(u,v)-expected)) < qptol failed.")
    end subroutine test_cross_product_rqp

    subroutine test_cross_product_csp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(sp) :: u(n),v(n),expected(n),diff(n)
        u = [cmplx(0,1,sp),cmplx(1,0,sp),cmplx(0,0,sp)]
        v = [cmplx(1,1,sp),cmplx(0,0,sp),cmplx(1,0,sp)]
        expected = [cmplx(1,0,sp),cmplx(0,-1,sp),cmplx(-1,-1,sp)]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < sptol), &
             "all(abs(cross_product(u,v)-expected)) < sptol failed.")
    end subroutine test_cross_product_csp

    subroutine test_cross_product_cdp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(dp) :: u(n),v(n),expected(n),diff(n)
        u = [cmplx(0,1,dp),cmplx(1,0,dp),cmplx(0,0,dp)]
        v = [cmplx(1,1,dp),cmplx(0,0,dp),cmplx(1,0,dp)]
        expected = [cmplx(1,0,dp),cmplx(0,-1,dp),cmplx(-1,-1,dp)]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < dptol), &
             "all(abs(cross_product(u,v)-expected)) < dptol failed.")
    end subroutine test_cross_product_cdp

    subroutine test_cross_product_cqp(error)
        !> Error handling
        type(error_type),allocatable,intent(out) :: error

        integer,parameter :: n = 3
        complex(qp) :: u(n),v(n),expected(n),diff(n)
        u = [cmplx(0,1,qp),cmplx(1,0,qp),cmplx(0,0,qp)]
        v = [cmplx(1,1,qp),cmplx(0,0,qp),cmplx(1,0,qp)]
        expected = [cmplx(1,0,qp),cmplx(0,-1,qp),cmplx(-1,-1,qp)]
        diff = expected - cross_product(u,v)
        call check(error,all(abs(diff) < qptol), &
             "all(abs(cross_product(u,v)-expected)) < qptol failed.")
    end subroutine test_cross_product_cqp

end module

program tester
    use,intrinsic :: iso_fortran_env,only:error_unit
    use testdrive,only:run_testsuite,new_testsuite,testsuite_type
    use test_linalg,only:collect_linalg
    implicit none
    integer :: stat,is
    type(testsuite_type),allocatable :: testsuites(:)
    character(len=*),parameter :: fmt = '("#", *(1x, a))'

    stat = 0

    testsuites = [ &
        new_testsuite("linalg",collect_linalg) &
        ]

    do is = 1,size(testsuites)
        write (error_unit,fmt) "Testing:",testsuites(is)%name
        call run_testsuite(testsuites(is)%collect,error_unit,stat)
    end do

    if (stat > 0) then
        write (error_unit,'(i0, 1x, a)') stat,"test(s) failed!"
        error stop
    end if
end program
