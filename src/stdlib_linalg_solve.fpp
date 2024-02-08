#:set REAL_KINDS    = ["sp", "dp", "qp"]
#:set REAL_INITIALS = ["s","d","q"]
#:set REAL_TYPES    = ["real({})".format(k) for k in REAL_KINDS]
#:set REAL_KINDS_TYPES = list(zip(REAL_KINDS, REAL_TYPES, REAL_INITIALS))

module stdlib_linalg_solve
     use stdlib_linalg_constants
     use stdlib_linalg_blas
     use stdlib_linalg_lapack
     use stdlib_linalg_state
     use iso_fortran_env,only:real32,real64,real128,int8,int16,int32,int64,stderr => error_unit
     implicit none(type,external)
     private

     !> Solve a linear system
     public :: solve


     ! Scipy: solve(a, b, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True, assume_a='gen', transposed=False)[source]#
     ! IMSL: lu_solve(a, b, transpose=False)


     interface solve
        #:for rk,rt,ri in REAL_KINDS_TYPES
        module procedure stdlib_linalg_${ri}$solve
        #:endfor
     end interface solve


     contains

     #:for rk,rt,ri in REAL_KINDS_TYPES
     ! Compute the solution to a real system of linear equations A * X = B
     function stdlib_linalg_${ri}$solve(a,b,err) result(x)
         real(${rk}$),                     intent(inout) :: a(:,:)
         real(${rk}$),                     intent(in)    :: b(:)
         type(linalg_state), optional, intent(out)   :: err
         real(${rk}$), allocatable, target :: x(:)

         ! Local variables
         type(linalg_state) :: err0
         integer(ilp) :: lda,n,ldb,nrhs,info
         integer(ilp), allocatable :: ipiv(:)
         real(${rk}$), pointer :: xmat(:,:)
         character(*), parameter :: this = 'solve'

         !> Problem sizes
         lda  = size(a,1,kind=ilp)
         n    = size(a,2,kind=ilp)
         ldb  = size(b,1,kind=ilp)
         nrhs = 1

         if (lda<1 .or. n<1 .or. ldb<1 .or. lda/=n .or. ldb/=n) then
            err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid sizes: a=[',lda,',',n,'], b=[',ldb,']')
            goto 1
         end if

         ! Pivot indices
         allocate(ipiv(n))

         ! Initialize solution with the rhs
         allocate(x(n),source=b)
         xmat(1:n,1:1) => x

         ! Solve system
         call gesv(n,nrhs,a,lda,ipiv,xmat,ldb,info)

         ! Process output
         select case (info)
            case (0)
                ! Success
            case (-1)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid rhs size n=',nrhs)
            case (-4)
                err0 = linalg_state(this,LINALG_VALUE_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (-7)
                err0 = linalg_state(this,LINALG_ERROR,'invalid matrix size a=[',lda,',',n,']')
            case (1:)
                err0 = linalg_state(this,LINALG_ERROR,'singular matrix')
            case default
                err0 = linalg_state(this,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

         ! Process output
         1 call linalg_error_handling(err0,err)

     end function stdlib_linalg_${ri}$solve


     #:endfor

end module stdlib_linalg_solve
