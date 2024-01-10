module stdlib_linalg_constants
     use iso_fortran_env, only: int32,int64
     use, intrinsic :: ieee_arithmetic, only: ieee_is_nan 
     implicit none(type,external)
     public



     integer, parameter :: sp = selected_real_kind(6)
     integer, parameter :: dp = selected_real_kind(15)
     integer, parameter :: lk = kind(.true.)



     ! 32-bit function prefixes 
     character,   parameter :: sprefix  = 'S' 
     character,   parameter :: cprefix  = 'C' 

     ! 32-bit real constants 
     real(sp),    parameter :: szero  =  0.0_sp
     real(sp),    parameter :: shalf  =  0.5_sp
     real(sp),    parameter :: sone   =  1.0_sp
     real(sp),    parameter :: stwo   =  2.0_sp
     real(sp),    parameter :: sthree =  3.0_sp
     real(sp),    parameter :: sfour  =  4.0_sp
     real(sp),    parameter :: seight =  8.0_sp
     real(sp),    parameter :: sten   = 10.0_sp

     ! 32-bit scaling constants 
     integer,     parameter :: smaxexp = maxexponent(szero) 
     integer,     parameter :: sminexp = minexponent(szero) 
     real(sp),    parameter :: sradix  = real(radix(szero),sp) 
     real(sp),    parameter :: sulp    = epsilon(szero) 
     real(sp),    parameter :: seps    = sulp*shalf 
     real(sp),    parameter :: ssafmin = sradix**max(sminexp-1,1-smaxexp) 
     real(sp),    parameter :: ssafmax = sone/ssafmin 
     real(sp),    parameter :: ssmlnum = ssafmin/sulp 
     real(sp),    parameter :: sbignum = ssafmax*sulp 
     real(sp),    parameter :: srtmin  = sqrt(ssmlnum) 
     real(sp),    parameter :: srtmax  = sqrt(sbignum) 

     ! 32-bit Blue's scaling constants 
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771 
     real(sp),    parameter :: stsml   = sradix**ceiling((sminexp-1)*shalf) 
     real(sp),    parameter :: stbig   = sradix**floor((smaxexp-digits(szero)+1)*shalf) 
     real(sp),    parameter :: sssml   = sradix**(-floor((sminexp-digits(szero))*shalf)) 
     real(sp),    parameter :: ssbig   = sradix**(-ceiling((smaxexp+digits(szero)-1)*shalf)) 

     ! 32-bit complex constants 
     complex(sp), parameter :: czero  = (0.0_sp,0.0_sp)
     complex(sp), parameter :: chalf  = (0.5_sp,0.0_sp)
     complex(sp), parameter :: cone   = (1.0_sp,0.0_sp)

     ! 64-bit function prefixes 
     character,   parameter :: dprefix  = 'D' 
     character,   parameter :: zprefix  = 'Z' 

     ! 64-bit real constants 
     real(dp),    parameter :: dzero  =  0.0_dp
     real(dp),    parameter :: dhalf  =  0.5_dp
     real(dp),    parameter :: done   =  1.0_dp
     real(dp),    parameter :: dtwo   =  2.0_dp
     real(dp),    parameter :: dthree =  3.0_dp
     real(dp),    parameter :: dfour  =  4.0_dp
     real(dp),    parameter :: deight =  8.0_dp
     real(dp),    parameter :: dten   = 10.0_dp

     ! 64-bit scaling constants 
     integer,     parameter :: dmaxexp = maxexponent(dzero) 
     integer,     parameter :: dminexp = minexponent(dzero) 
     real(dp),    parameter :: dradix  = real(radix(dzero),dp) 
     real(dp),    parameter :: dulp    = epsilon(dzero) 
     real(dp),    parameter :: deps    = dulp*dhalf 
     real(dp),    parameter :: dsafmin = dradix**max(dminexp-1,1-dmaxexp) 
     real(dp),    parameter :: dsafmax = done/dsafmin 
     real(dp),    parameter :: dsmlnum = dsafmin/dulp 
     real(dp),    parameter :: dbignum = dsafmax*dulp 
     real(dp),    parameter :: drtmin  = sqrt(dsmlnum) 
     real(dp),    parameter :: drtmax  = sqrt(dbignum) 

     ! 64-bit Blue's scaling constants 
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771 
     real(dp),    parameter :: dtsml   = dradix**ceiling((dminexp-1)*dhalf) 
     real(dp),    parameter :: dtbig   = dradix**floor((dmaxexp-digits(dzero)+1)*dhalf) 
     real(dp),    parameter :: dssml   = dradix**(-floor((dminexp-digits(dzero))*dhalf)) 
     real(dp),    parameter :: dsbig   = dradix**(-ceiling((dmaxexp+digits(dzero)-1)*dhalf)) 

     ! 64-bit complex constants 
     complex(dp), parameter :: zzero  = (0.0_dp,0.0_dp)
     complex(dp), parameter :: zhalf  = (0.5_dp,0.0_dp)
     complex(dp), parameter :: zone   = (1.0_dp,0.0_dp)



contains







end module stdlib_linalg_constants
