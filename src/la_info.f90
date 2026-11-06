!> @brief Decode the `info` exit code of the LAPACK routines into a `la_state`.
!>
!> One `elemental` handler per LAPACK routine, `handle_<routine>`, shared by every procedure
!> that calls it. Each handler takes the caller's name first, which becomes the `where_at` of
!> the returned state, then the exit code, then the routine arguments its messages need.
!>
!> One classification holds for every routine: an invalid dimension of the user's data is a
!> `LINALG_VALUE_ERROR`; an invalid argument the library chose (a job character, a triangle
!> selection, a workspace size) is a `LINALG_INTERNAL_ERROR`; a failure that `info > 0`
!> reports is a `LINALG_ERROR`; any code that LAPACK does not document is a
!> `LINALG_INTERNAL_ERROR`.
module la_info
     use la_constants,only:ilp
     use la_state_type,only:la_state,LINALG_ERROR,LINALG_INTERNAL_ERROR,LINALG_VALUE_ERROR
     implicit none(type,external)
     private

     public :: handle_gesv
     public :: handle_getrf
     public :: handle_getri
     public :: handle_gecon
     public :: handle_posv
     public :: handle_potrs
     public :: handle_potrf
     public :: handle_geqrf
     public :: handle_geqp3
     public :: handle_orgqr
     public :: handle_gees
     public :: handle_geev
     public :: handle_ggev
     public :: handle_heev
     public :: handle_gesdd
     public :: handle_gelsd
     public :: handle_gglse
     public :: handle_ggglm

     contains

     !> Process the exit code of gesv(n,nrhs,a,lda,ipiv,b,ldb,info)
     elemental subroutine handle_gesv(where_at,info,n,nrhs,lda,ldb,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gesv exit code
         integer(ilp), intent(in) :: info
         !> gesv arguments
         integer(ilp), intent(in) :: n,nrhs,lda,ldb
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid problem size n=',n)
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid rhs size nrhs=',nrhs)
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size a=',[lda,n])
            case (-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid rhs size b=',[ldb,nrhs])
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'singular matrix')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gesv

     !> Process the exit code of getrf(m,n,a,lda,ipiv,info)
     elemental subroutine handle_getrf(where_at,info,m,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> getrf exit code
         integer(ilp), intent(in) :: info
         !> getrf arguments
         integer(ilp), intent(in) :: m,n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1,-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size a=[',m,',',n,']')
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',m)
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'singular matrix')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_getrf

     !> Process the exit code of getri(n,a,lda,ipiv,work,lwork,info)
     elemental subroutine handle_getri(where_at,info,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> getri exit code
         integer(ilp), intent(in) :: info
         !> getri arguments
         integer(ilp), intent(in) :: n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-3)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',n)
            case (-6)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size')
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'singular matrix')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_getri

     !> Process the exit code of gecon(norm,n,a,lda,anorm,rcond,work,iwork,info)
     elemental subroutine handle_gecon(where_at,info,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gecon exit code
         integer(ilp), intent(in) :: info
         !> gecon arguments
         integer(ilp), intent(in) :: n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid norm type')
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',n)
            case (-5)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid norm of a: negative, NaN or infinite')
            case (1)
                err = la_state(where_at,LINALG_ERROR,'the reciprocal condition number is NaN or infinite, ', &
                                                     'or the inverse has zero norm')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gecon

     !> Process the exit code of posv(uplo,n,nrhs,a,lda,b,ldb,info)
     elemental subroutine handle_posv(where_at,info,triangle,n,nrhs,lda,ldb,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> posv exit code
         integer(ilp), intent(in) :: info
         !> posv arguments
         character, intent(in) :: triangle
         integer(ilp), intent(in) :: n,nrhs,lda,ldb
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid triangle selection: ', &
                                    triangle,'. should be U/L')
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-3)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid rhs size nrhs=',nrhs)
            case (-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',n)
            case (-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid ldb=',ldb,': should be >=',n)
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'matrix is not positive definite: ', &
                                    'leading minor of order',info,'is not positive definite')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_posv

     !> Process the exit code of potrs(uplo,n,nrhs,a,lda,b,ldb,info)
     elemental subroutine handle_potrs(where_at,info,triangle,n,nrhs,lda,ldb,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> potrs exit code
         integer(ilp), intent(in) :: info
         !> potrs arguments
         character, intent(in) :: triangle
         integer(ilp), intent(in) :: n,nrhs,lda,ldb
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid triangle selection: ', &
                                    triangle,'. should be U/L')
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-3)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid rhs size nrhs=',nrhs)
            case (-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': should be >=',n)
            case (-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid ldb=',ldb,': should be >=',n)
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_potrs

     !> Process the exit code of potrf(uplo,n,a,lda,info)
     elemental subroutine handle_potrf(where_at,info,triangle,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> potrf exit code
         integer(ilp), intent(in) :: info
         !> potrf arguments
         character, intent(in) :: triangle
         integer(ilp), intent(in) :: n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid triangle selection: ', &
                                   triangle,'. should be U/L')
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid lda=',lda,': is < n = ',n)
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'cannot complete factorization:',info, &
                                   '-th order leading minor is not positive definite')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_potrf

     !> Process the exit code of geqrf(m,n,a,lda,tau,work,lwork,info)
     elemental subroutine handle_geqrf(where_at,info,m,n,lwork,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> geqrf exit code
         integer(ilp), intent(in) :: info
         !> geqrf arguments
         integer(ilp), intent(in) :: m,n,lwork
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size m=',m)
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size a=',[m,n])
            case (-7)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid input for lwork=',lwork)
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_geqrf

     !> Process the exit code of geqp3(m,n,a,lda,jpvt,tau,work,lwork,[rwork,]info)
     elemental subroutine handle_geqp3(where_at,info,m,n,lwork,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> geqp3 exit code
         integer(ilp), intent(in) :: info
         !> geqp3 arguments
         integer(ilp), intent(in) :: m,n,lwork
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size m=',m)
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size a=',[m,n])
            case (-8)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid input for lwork=',lwork)
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_geqp3

     !> Process the exit code of orgqr/ungqr(m,n,k,a,lda,tau,work,lwork,info)
     elemental subroutine handle_orgqr(where_at,info,m,n,k,lwork,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> orgqr exit code
         integer(ilp), intent(in) :: info
         !> orgqr arguments
         integer(ilp), intent(in) :: m,n,k,lwork
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size m=',m)
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size n=',n)
            case (-3)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid k=min(m,n)=',k)
            case (-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size a=',[m,n])
            case (-8)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'invalid input for lwork=',lwork)
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_orgqr

     !> Process the exit code of gees: real gees(jobvs,sort,select,n,a,lda,sdim,wr,wi,vs,ldvs,work,
     !> lwork,bwork,info), complex gees(jobvs,sort,select,n,a,lda,sdim,w,vs,ldvs,work,lwork,rwork,
     !> bwork,info). The complex form has one eigenvalue array less, so ldvs and lwork come one
     !> position earlier: never on a code where the real form has a scalar argument
     elemental subroutine handle_gees(where_at,info,n,lda,ldvs,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gees exit code
         integer(ilp), intent(in) :: info
         !> gees arguments
         integer(ilp), intent(in) :: n,lda,ldvs
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid Schur vector task request')
            case (-2)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid sorting task request')
            case (-4,-6)
                err = la_state(where_at,LINALG_VALUE_ERROR,'Invalid/non-square input matrix size:',[lda,n])
            case (-10,-11)
                err = la_state(where_at,LINALG_VALUE_ERROR,'Schur vector matrix has insufficient size',[ldvs,n])
            case (-12,-13)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Insufficient working storage size')
            case (1:)
                if (info<n) then
                    err = la_state(where_at,LINALG_ERROR,'Convergence failure; converged range is',[info,n])
                elseif (info==n) then
                    err = la_state(where_at,LINALG_ERROR,'Convergence failure: no converged eigenvalues')
                elseif (info==n+1) then
                    err = la_state(where_at,LINALG_ERROR,'Ill-conditioned problem: could not sort eigenvalues')
                elseif (info==n+2) then
                    err = la_state(where_at,LINALG_ERROR,'Some selected eigenvalues lost property due to sorting')
                else
                    err = la_state(where_at,LINALG_INTERNAL_ERROR,'GEES catastrophic error: info=',info)
                end if
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'GEES catastrophic error: info=',info)
         end select

     end subroutine handle_gees

     !> Process the exit code of geev: real geev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,
     !> info), complex geev(jobvl,jobvr,n,a,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info). The complex
     !> form has one eigenvalue array less, so ldvl, ldvr and lwork come one position earlier: never
     !> on a code where the real form has a scalar argument
     elemental subroutine handle_geev(where_at,info,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> geev exit code
         integer(ilp), intent(in) :: info
         !> geev arguments
         integer(ilp), intent(in) :: n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid task ID: left eigenvectors.')
            case (-2)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid task ID: right eigenvectors.')
            case (-3,-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size: a=',[lda,n])
            case (-8,-9)
                err = la_state(where_at,LINALG_VALUE_ERROR,'insufficient left vector matrix size.')
            case (-10,-11)
                err = la_state(where_at,LINALG_VALUE_ERROR,'insufficient right vector matrix size.')
            case (-12,-13)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Insufficient work array size.')
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'Eigenvalue computation did not converge.')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Unknown error returned by geev.')
         end select

     end subroutine handle_geev

     !> Process the exit code of ggev: real ggev(jobvl,jobvr,n,a,lda,b,ldb,alphar,alphai,beta,vl,ldvl,
     !> vr,ldvr,work,lwork,info), complex ggev(jobvl,jobvr,n,a,lda,b,ldb,alpha,beta,vl,ldvl,vr,ldvr,
     !> work,lwork,rwork,info). The complex form has one eigenvalue array less, so ldvl, ldvr and
     !> lwork come one position earlier: never on a code where the real form has a scalar argument
     elemental subroutine handle_ggev(where_at,info,n,lda,ldb,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> ggev exit code
         integer(ilp), intent(in) :: info
         !> ggev arguments
         integer(ilp), intent(in) :: n,lda,ldb
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid task ID: left eigenvectors.')
            case (-2)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid task ID: right eigenvectors.')
            case (-3,-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size: a=',[lda,n])
            case (-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size: b=',[ldb,n])
            case (-11,-12)
                err = la_state(where_at,LINALG_VALUE_ERROR,'insufficient left vector matrix size.')
            case (-13,-14)
                err = la_state(where_at,LINALG_VALUE_ERROR,'insufficient right vector matrix size.')
            case (-15,-16)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Insufficient work array size.')
            case (1:)
                if (info<=n) then
                    err = la_state(where_at,LINALG_ERROR,'Eigenvalue computation did not converge.')
                elseif (info==n+1) then
                    err = la_state(where_at,LINALG_ERROR,'Eigenvalue computation failed.')
                elseif (info==n+2) then
                    err = la_state(where_at,LINALG_ERROR,'Eigenvector computation failed.')
                else
                    err = la_state(where_at,LINALG_INTERNAL_ERROR,'Unknown error returned by ggev.')
                end if
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Unknown error returned by ggev.')
         end select

     end subroutine handle_ggev

     !> Process the exit code of syev(jobz,uplo,n,a,lda,w,work,lwork,info) and
     !> heev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
     elemental subroutine handle_heev(where_at,info,n,lda,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> syev/heev exit code
         integer(ilp), intent(in) :: info
         !> syev/heev arguments
         integer(ilp), intent(in) :: n,lda
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid eigenvector request.')
            case (-2)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid triangular section request.')
            case (-3,-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size: a=',[lda,n])
            case (-8)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size.')
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'Eigenvalue computation did not converge.')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Unknown error returned by syev/heev.')
         end select

     end subroutine handle_heev

     !> Process the exit code of gesdd(jobz,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,[rwork,]iwork,info)
     elemental subroutine handle_gesdd(where_at,info,m,n,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gesdd exit code
         integer(ilp), intent(in) :: info
         !> gesdd arguments
         integer(ilp), intent(in) :: m,n
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Invalid task ID on input to GESDD.')
            case (-2,-3,-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix size: a=[',m,',',n,']')
            case (-4)
                err = la_state(where_at,LINALG_VALUE_ERROR,'A contains invalid/NaN values.')
            case (-8)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix U size, with a=[',m,',',n,']')
            case (-10)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid matrix V size, with a=[',m,',',n,']')
            case (-12)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size')
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'SVD computation did not converge.')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'Unknown error returned by GESDD.')
         end select

     end subroutine handle_gesdd

     !> Process the exit code of gelsd(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,[rwork,]iwork,info)
     elemental subroutine handle_gelsd(where_at,info,n,nrhs,lda,ldb,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gelsd exit code
         integer(ilp), intent(in) :: info
         !> gelsd arguments
         integer(ilp), intent(in) :: n,nrhs,lda,ldb
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (-1,-2,-3,-5,-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid problem size a=',[lda,n],', b=',[ldb,nrhs])
            case (-12)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size')
            case (1:)
                err = la_state(where_at,LINALG_ERROR,'SVD did not converge.')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gelsd

     !> Process the exit code of gglse(m,n,p,a,lda,b,ldb,c,d,x,work,lwork,info), that minimizes
     !> ||c - a*x|| subject to b*x = d. The messages name the matrices as constrained_lstsq does:
     !> a(m,n) and the constraint c(p,n)
     elemental subroutine handle_gglse(where_at,info,m,n,p,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> gglse exit code
         integer(ilp), intent(in) :: info
         !> gglse arguments
         integer(ilp), intent(in) :: m,n,p
         !> State return flag
         type(la_state), intent(out) :: err

         select case (info)
            case (0)
                ! Success
            case (1)
                err = la_state(where_at,LINALG_ERROR,'rank(c) < p: the least-squares solution cannot be computed')
            case (2)
                err = la_state(where_at,LINALG_ERROR,'rank([a; c]) < n: the least-squares solution cannot be computed')
            case (-1,-5)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid number of rows for a, m=',m)
            case (-2)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid number of columns for a and c, n=',n)
            case (-3,-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid number of rows for c, p=',p)
            case (-12)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_gglse

     !> Process the exit code of ggglm(n,m,p,a,lda,b,ldb,d,x,y,work,lwork,info), that minimizes
     !> ||y|| subject to d = a*x + b*y. The messages name the matrices as generalized_lstsq does:
     !> a(m,n), with m rows, and the square root w of the covariance in place of b
     elemental subroutine handle_ggglm(where_at,info,m,n,err)
         !> Name of the calling procedure
         character(*), intent(in) :: where_at
         !> ggglm exit code
         integer(ilp), intent(in) :: info
         !> Rows and columns of a, ggglm's n and m
         integer(ilp), intent(in) :: m,n
         !> State return flag
         type(la_state), intent(out) :: err

         ! ggglm returns 1 when the solve with T22 fails (rank([a w]) < m), and 2 when the solve with
         ! R11 fails (rank(a) < n)
         select case (info)
            case (0)
                ! Success
            case (1)
                err = la_state(where_at,LINALG_ERROR,'the square root of w is singular')
            case (2)
                err = la_state(where_at,LINALG_ERROR,'a does not have full column rank')
            case (-1,-2,-3,-5,-7)
                err = la_state(where_at,LINALG_VALUE_ERROR,'invalid problem size a(m,n) =',[m,n])
            case (-12)
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'insufficient workspace size')
            case default
                err = la_state(where_at,LINALG_INTERNAL_ERROR,'catastrophic error')
         end select

     end subroutine handle_ggglm

end module la_info
