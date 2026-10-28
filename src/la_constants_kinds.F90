!> 32-bit BLAS/LAPACK constants
module la_constants_sp
     use la_constants
     implicit none(type,external)
     private

     public :: negone,zero,half,one,two,three,four,eight,ten
     public :: czero,chalf,cone,cnegone
     public :: maxexp,minexp,rradix,ulp,eps,safmin,safmax,smlnum,bignum,rtmin,rtmax
     public :: tsml,tbig,ssml,sbig

     ! 32-bit real constants
     real(sp),parameter :: negone = -1.00_sp
     real(sp),parameter :: zero = 0.00_sp
     real(sp),parameter :: half = 0.50_sp
     real(sp),parameter :: one = 1.00_sp
     real(sp),parameter :: two = 2.00_sp
     real(sp),parameter :: three = 3.00_sp
     real(sp),parameter :: four = 4.00_sp
     real(sp),parameter :: eight = 8.00_sp
     real(sp),parameter :: ten = 10.00_sp

     ! 32-bit complex constants
     complex(sp),parameter :: czero = (0.0_sp,0.0_sp)
     complex(sp),parameter :: chalf = (0.5_sp,0.0_sp)
     complex(sp),parameter :: cone = (1.0_sp,0.0_sp)
     complex(sp),parameter :: cnegone = (-1.0_sp,0.0_sp)

     ! 32-bit scaling constants
     integer,parameter :: maxexp = maxexponent(zero)
     integer,parameter :: minexp = minexponent(zero)
     real(sp),parameter :: rradix = real(radix(zero),sp)
     real(sp),parameter :: ulp = epsilon(zero)
     real(sp),parameter :: eps = ulp*half
     real(sp),parameter :: safmin = rradix**max(minexp - 1,1 - maxexp)
     real(sp),parameter :: safmax = one/safmin
     real(sp),parameter :: smlnum = safmin/ulp
     real(sp),parameter :: bignum = safmax*ulp
     real(sp),parameter :: rtmin = sqrt(smlnum)
     real(sp),parameter :: rtmax = sqrt(bignum)

     ! 32-bit Blue's scaling constants
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771
     real(sp),parameter :: tsml = rradix**ceiling((minexp - 1)*half)
     real(sp),parameter :: tbig = rradix**floor((maxexp - digits(zero) + 1)*half)
     real(sp),parameter :: ssml = rradix**(-floor((minexp - digits(zero))*half))
     real(sp),parameter :: sbig = rradix**(-ceiling((maxexp + digits(zero) - 1)*half))

end module la_constants_sp

!> 64-bit BLAS/LAPACK constants
module la_constants_dp
     use la_constants
     implicit none(type,external)
     private

     public :: negone,zero,half,one,two,three,four,eight,ten
     public :: czero,chalf,cone,cnegone
     public :: maxexp,minexp,rradix,ulp,eps,safmin,safmax,smlnum,bignum,rtmin,rtmax
     public :: tsml,tbig,ssml,sbig

     ! 64-bit real constants
     real(dp),parameter :: negone = -1.00_dp
     real(dp),parameter :: zero = 0.00_dp
     real(dp),parameter :: half = 0.50_dp
     real(dp),parameter :: one = 1.00_dp
     real(dp),parameter :: two = 2.00_dp
     real(dp),parameter :: three = 3.00_dp
     real(dp),parameter :: four = 4.00_dp
     real(dp),parameter :: eight = 8.00_dp
     real(dp),parameter :: ten = 10.00_dp

     ! 64-bit complex constants
     complex(dp),parameter :: czero = (0.0_dp,0.0_dp)
     complex(dp),parameter :: chalf = (0.5_dp,0.0_dp)
     complex(dp),parameter :: cone = (1.0_dp,0.0_dp)
     complex(dp),parameter :: cnegone = (-1.0_dp,0.0_dp)

     ! 64-bit scaling constants
     integer,parameter :: maxexp = maxexponent(zero)
     integer,parameter :: minexp = minexponent(zero)
     real(dp),parameter :: rradix = real(radix(zero),dp)
     real(dp),parameter :: ulp = epsilon(zero)
     real(dp),parameter :: eps = ulp*half
     real(dp),parameter :: safmin = rradix**max(minexp - 1,1 - maxexp)
     real(dp),parameter :: safmax = one/safmin
     real(dp),parameter :: smlnum = safmin/ulp
     real(dp),parameter :: bignum = safmax*ulp
     real(dp),parameter :: rtmin = sqrt(smlnum)
     real(dp),parameter :: rtmax = sqrt(bignum)

     ! 64-bit Blue's scaling constants
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771
     real(dp),parameter :: tsml = rradix**ceiling((minexp - 1)*half)
     real(dp),parameter :: tbig = rradix**floor((maxexp - digits(zero) + 1)*half)
     real(dp),parameter :: ssml = rradix**(-floor((minexp - digits(zero))*half))
     real(dp),parameter :: sbig = rradix**(-ceiling((maxexp + digits(zero) - 1)*half))

end module la_constants_dp

#ifdef LA_WITH_XDP
!> 80-bit BLAS/LAPACK constants
module la_constants_xdp
     use la_constants
     implicit none(type,external)
     private

     public :: negone,zero,half,one,two,three,four,eight,ten
     public :: czero,chalf,cone,cnegone
     public :: maxexp,minexp,rradix,ulp,eps,safmin,safmax,smlnum,bignum,rtmin,rtmax
     public :: tsml,tbig,ssml,sbig

     ! 80-bit real constants
     real(xdp),parameter :: negone = -1.00_xdp
     real(xdp),parameter :: zero = 0.00_xdp
     real(xdp),parameter :: half = 0.50_xdp
     real(xdp),parameter :: one = 1.00_xdp
     real(xdp),parameter :: two = 2.00_xdp
     real(xdp),parameter :: three = 3.00_xdp
     real(xdp),parameter :: four = 4.00_xdp
     real(xdp),parameter :: eight = 8.00_xdp
     real(xdp),parameter :: ten = 10.00_xdp

     ! 80-bit complex constants
     complex(xdp),parameter :: czero = (0.0_xdp,0.0_xdp)
     complex(xdp),parameter :: chalf = (0.5_xdp,0.0_xdp)
     complex(xdp),parameter :: cone = (1.0_xdp,0.0_xdp)
     complex(xdp),parameter :: cnegone = (-1.0_xdp,0.0_xdp)

     ! 80-bit scaling constants
     integer,parameter :: maxexp = maxexponent(zero)
     integer,parameter :: minexp = minexponent(zero)
     real(xdp),parameter :: rradix = real(radix(zero),xdp)
     real(xdp),parameter :: ulp = epsilon(zero)
     real(xdp),parameter :: eps = ulp*half
     real(xdp),parameter :: safmin = rradix**max(minexp - 1,1 - maxexp)
     real(xdp),parameter :: safmax = one/safmin
     real(xdp),parameter :: smlnum = safmin/ulp
     real(xdp),parameter :: bignum = safmax*ulp
     real(xdp),parameter :: rtmin = sqrt(smlnum)
     real(xdp),parameter :: rtmax = sqrt(bignum)

     ! 80-bit Blue's scaling constants
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771
     real(xdp),parameter :: tsml = rradix**ceiling((minexp - 1)*half)
     real(xdp),parameter :: tbig = rradix**floor((maxexp - digits(zero) + 1)*half)
     real(xdp),parameter :: ssml = rradix**(-floor((minexp - digits(zero))*half))
     real(xdp),parameter :: sbig = rradix**(-ceiling((maxexp + digits(zero) - 1)*half))

end module la_constants_xdp
#endif

#ifdef LA_WITH_QP
!> 128-bit BLAS/LAPACK constants
module la_constants_qp
     use la_constants
     implicit none(type,external)
     private

     public :: negone,zero,half,one,two,three,four,eight,ten
     public :: czero,chalf,cone,cnegone
     public :: maxexp,minexp,rradix,ulp,eps,safmin,safmax,smlnum,bignum,rtmin,rtmax
     public :: tsml,tbig,ssml,sbig

     ! 128-bit real constants
     real(qp),parameter :: negone = -1.00_qp
     real(qp),parameter :: zero = 0.00_qp
     real(qp),parameter :: half = 0.50_qp
     real(qp),parameter :: one = 1.00_qp
     real(qp),parameter :: two = 2.00_qp
     real(qp),parameter :: three = 3.00_qp
     real(qp),parameter :: four = 4.00_qp
     real(qp),parameter :: eight = 8.00_qp
     real(qp),parameter :: ten = 10.00_qp

     ! 128-bit complex constants
     complex(qp),parameter :: czero = (0.0_qp,0.0_qp)
     complex(qp),parameter :: chalf = (0.5_qp,0.0_qp)
     complex(qp),parameter :: cone = (1.0_qp,0.0_qp)
     complex(qp),parameter :: cnegone = (-1.0_qp,0.0_qp)

     ! 128-bit scaling constants
     integer,parameter :: maxexp = maxexponent(zero)
     integer,parameter :: minexp = minexponent(zero)
     real(qp),parameter :: rradix = real(radix(zero),qp)
     real(qp),parameter :: ulp = epsilon(zero)
     real(qp),parameter :: eps = ulp*half
     real(qp),parameter :: safmin = rradix**max(minexp - 1,1 - maxexp)
     real(qp),parameter :: safmax = one/safmin
     real(qp),parameter :: smlnum = safmin/ulp
     real(qp),parameter :: bignum = safmax*ulp
     real(qp),parameter :: rtmin = sqrt(smlnum)
     real(qp),parameter :: rtmax = sqrt(bignum)

     ! 128-bit Blue's scaling constants
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771
     real(qp),parameter :: tsml = rradix**ceiling((minexp - 1)*half)
     real(qp),parameter :: tbig = rradix**floor((maxexp - digits(zero) + 1)*half)
     real(qp),parameter :: ssml = rradix**(-floor((minexp - digits(zero))*half))
     real(qp),parameter :: sbig = rradix**(-ceiling((maxexp + digits(zero) - 1)*half))

end module la_constants_qp
#endif

