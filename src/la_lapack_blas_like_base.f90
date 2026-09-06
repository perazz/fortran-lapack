!> BLAS-like base: copy, precision conversion, random and packed storage
module la_lapack_blas_like_base
     use la_constants
     use la_blas_aux
     use la_lapack_auxiliary
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slacpy
     public :: la_slaruv
     public :: la_slaset
     public :: la_stfttp
     public :: la_stfttr
     public :: la_stpttf
     public :: la_stpttr
     public :: la_strttf
     public :: la_strttp
     public :: la_slag2d
     public :: la_slarnv
     public :: la_dlacpy
     public :: la_dlag2s
     public :: la_dlaruv
     public :: la_dlaset
     public :: la_dlat2s
     public :: la_dtfttp
     public :: la_dtfttr
     public :: la_dtpttf
     public :: la_dtpttr
     public :: la_dtrttf
     public :: la_dtrttp
     public :: la_dlag2q
     public :: la_dlarnv
     public :: la_qlacpy
     public :: la_qlag2d
     public :: la_qlaruv
     public :: la_qlaset
     public :: la_qlat2d
     public :: la_qtfttp
     public :: la_qtfttr
     public :: la_qtpttf
     public :: la_qtpttr
     public :: la_qtrttf
     public :: la_qtrttp
     public :: la_qlarnv
     public :: la_clag2z
     public :: la_clacp2
     public :: la_clacpy
     public :: la_clarnv
     public :: la_claset
     public :: la_ctfttp
     public :: la_ctfttr
     public :: la_ctpttf
     public :: la_ctpttr
     public :: la_ctrttf
     public :: la_ctrttp
     public :: la_zlag2w
     public :: la_zlacp2
     public :: la_zlacpy
     public :: la_zlag2c
     public :: la_zlarnv
     public :: la_zlaset
     public :: la_zlat2c
     public :: la_ztfttp
     public :: la_ztfttr
     public :: la_ztpttf
     public :: la_ztpttr
     public :: la_ztrttf
     public :: la_ztrttp
     public :: la_wlacp2
     public :: la_wlacpy
     public :: la_wlag2z
     public :: la_wlarnv
     public :: la_wlaset
     public :: la_wlat2z
     public :: la_wtfttp
     public :: la_wtfttr
     public :: la_wtpttf
     public :: la_wtpttr
     public :: la_wtrttf
     public :: la_wtrttp

     contains

     !> SLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_slacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_slacpy
     !> DLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_dlacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_dlacpy
     !> QLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_qlacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_qlacpy

     !> DLAG2S: converts a DOUBLE PRECISION matrix, SA, to a SINGLE
     !> PRECISION matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> DLAG2S checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_dlag2s(m,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           real(sp),intent(out) :: sa(ldsa,*)
           real(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: rmax
           ! Executable Statements
           rmax = la_slamch('O')
           do j = 1,n
              do i = 1,m
                 if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                    info = 1
                    go to 30
                 end if
                 sa(i,j) = a(i,j)
              end do
           end do
           info = 0
           30 continue
           return
     end subroutine la_dlag2s
     !> QLAG2D: converts a QUAD PRECISION matrix, SA, to a SINGLE
     !> PRECISION matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> QLAG2D checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_qlag2d(m,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           real(dp),intent(out) :: sa(ldsa,*)
           real(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: rmax
           ! Executable Statements
           rmax = la_dlamch('O')
           do j = 1,n
              do i = 1,m
                 if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                    info = 1
                    go to 30
                 end if
                 sa(i,j) = a(i,j)
              end do
           end do
           info = 0
           30 continue
           return
     end subroutine la_qlag2d

     !> SLARUV: returns a vector of n random real numbers from a uniform (0,1)
     !> distribution (n <= 128).
     !> This is an auxiliary routine called by SLARNV and CLARNV.

     pure subroutine la_slaruv(iseed,n,x)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(sp),intent(out) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           integer(ilp),parameter :: ipw2 = 4096
           real(sp),parameter :: r = one/ipw2

           ! Local Scalars
           integer(ilp) :: i,i1,i2,i3,i4,it1,it2,it3,it4,j
           ! Local Arrays
           integer(ilp) :: mm(lv,4)
           ! Intrinsic Functions
           intrinsic :: min,mod,real
           ! Data Statements
           mm(1,1:4) = [494,322,2508,2549]
           mm(2,1:4) = [2637,789,3754,1145]
           mm(3,1:4) = [255,1440,1766,2253]
           mm(4,1:4) = [2008,752,3572,305]
           mm(5,1:4) = [1253,2859,2893,3301]
           mm(6,1:4) = [3344,123,307,1065]
           mm(7,1:4) = [4084,1848,1297,3133]
           mm(8,1:4) = [1739,643,3966,2913]
           mm(9,1:4) = [3143,2405,758,3285]
           mm(10,1:4) = [3468,2638,2598,1241]
           mm(11,1:4) = [688,2344,3406,1197]
           mm(12,1:4) = [1657,46,2922,3729]
           mm(13,1:4) = [1238,3814,1038,2501]
           mm(14,1:4) = [3166,913,2934,1673]
           mm(15,1:4) = [1292,3649,2091,541]
           mm(16,1:4) = [3422,339,2451,2753]
           mm(17,1:4) = [1270,3808,1580,949]
           mm(18,1:4) = [2016,822,1958,2361]
           mm(19,1:4) = [154,2832,2055,1165]
           mm(20,1:4) = [2862,3078,1507,4081]
           mm(21,1:4) = [697,3633,1078,2725]
           mm(22,1:4) = [1706,2970,3273,3305]
           mm(23,1:4) = [491,637,17,3069]
           mm(24,1:4) = [931,2249,854,3617]
           mm(25,1:4) = [1444,2081,2916,3733]
           mm(26,1:4) = [444,4019,3971,409]
           mm(27,1:4) = [3577,1478,2889,2157]
           mm(28,1:4) = [3944,242,3831,1361]
           mm(29,1:4) = [2184,481,2621,3973]
           mm(30,1:4) = [1661,2075,1541,1865]
           mm(31,1:4) = [3482,4058,893,2525]
           mm(32,1:4) = [657,622,736,1409]
           mm(33,1:4) = [3023,3376,3992,3445]
           mm(34,1:4) = [3618,812,787,3577]
           mm(35,1:4) = [1267,234,2125,77]
           mm(36,1:4) = [1828,641,2364,3761]
           mm(37,1:4) = [164,4005,2460,2149]
           mm(38,1:4) = [3798,1122,257,1449]
           mm(39,1:4) = [3087,3135,1574,3005]
           mm(40,1:4) = [2400,2640,3912,225]
           mm(41,1:4) = [2870,2302,1216,85]
           mm(42,1:4) = [3876,40,3248,3673]
           mm(43,1:4) = [1905,1832,3401,3117]
           mm(44,1:4) = [1593,2247,2124,3089]
           mm(45,1:4) = [1797,2034,2762,1349]
           mm(46,1:4) = [1234,2637,149,2057]
           mm(47,1:4) = [3460,1287,2245,413]
           mm(48,1:4) = [328,1691,166,65]
           mm(49,1:4) = [2861,496,466,1845]
           mm(50,1:4) = [1950,1597,4018,697]
           mm(51,1:4) = [617,2394,1399,3085]
           mm(52,1:4) = [2070,2584,190,3441]
           mm(53,1:4) = [3331,1843,2879,1573]
           mm(54,1:4) = [769,336,153,3689]
           mm(55,1:4) = [1558,1472,2320,2941]
           mm(56,1:4) = [2412,2407,18,929]
           mm(57,1:4) = [2800,433,712,533]
           mm(58,1:4) = [189,2096,2159,2841]
           mm(59,1:4) = [287,1761,2318,4077]
           mm(60,1:4) = [2045,2810,2091,721]
           mm(61,1:4) = [1227,566,3443,2821]
           mm(62,1:4) = [2838,442,1510,2249]
           mm(63,1:4) = [209,41,449,2397]
           mm(64,1:4) = [2770,1238,1956,2817]
           mm(65,1:4) = [3654,1086,2201,245]
           mm(66,1:4) = [3993,603,3137,1913]
           mm(67,1:4) = [192,840,3399,1997]
           mm(68,1:4) = [2253,3168,1321,3121]
           mm(69,1:4) = [3491,1499,2271,997]
           mm(70,1:4) = [2889,1084,3667,1833]
           mm(71,1:4) = [2857,3438,2703,2877]
           mm(72,1:4) = [2094,2408,629,1633]
           mm(73,1:4) = [1818,1589,2365,981]
           mm(74,1:4) = [688,2391,2431,2009]
           mm(75,1:4) = [1407,288,1113,941]
           mm(76,1:4) = [634,26,3922,2449]
           mm(77,1:4) = [3231,512,2554,197]
           mm(78,1:4) = [815,1456,184,2441]
           mm(79,1:4) = [3524,171,2099,285]
           mm(80,1:4) = [1914,1677,3228,1473]
           mm(81,1:4) = [516,2657,4012,2741]
           mm(82,1:4) = [164,2270,1921,3129]
           mm(83,1:4) = [303,2587,3452,909]
           mm(84,1:4) = [2144,2961,3901,2801]
           mm(85,1:4) = [3480,1970,572,421]
           mm(86,1:4) = [119,1817,3309,4073]
           mm(87,1:4) = [3357,676,3171,2813]
           mm(88,1:4) = [837,1410,817,2337]
           mm(89,1:4) = [2826,3723,3039,1429]
           mm(90,1:4) = [2332,2803,1696,1177]
           mm(91,1:4) = [2089,3185,1256,1901]
           mm(92,1:4) = [3780,184,3715,81]
           mm(93,1:4) = [1700,663,2077,1669]
           mm(94,1:4) = [3712,499,3019,2633]
           mm(95,1:4) = [150,3784,1497,2269]
           mm(96,1:4) = [2000,1631,1101,129]
           mm(97,1:4) = [3375,1925,717,1141]
           mm(98,1:4) = [1621,3912,51,249]
           mm(99,1:4) = [3090,1398,981,3917]
           mm(100,1:4) = [3765,1349,1978,2481]
           mm(101,1:4) = [1149,1441,1813,3941]
           mm(102,1:4) = [3146,2224,3881,2217]
           mm(103,1:4) = [33,2411,76,2749]
           mm(104,1:4) = [3082,1907,3846,3041]
           mm(105,1:4) = [2741,3192,3694,1877]
           mm(106,1:4) = [359,2786,1682,345]
           mm(107,1:4) = [3316,382,124,2861]
           mm(108,1:4) = [1749,37,1660,1809]
           mm(109,1:4) = [185,759,3997,3141]
           mm(110,1:4) = [2784,2948,479,2825]
           mm(111,1:4) = [2202,1862,1141,157]
           mm(112,1:4) = [2199,3802,886,2881]
           mm(113,1:4) = [1364,2423,3514,3637]
           mm(114,1:4) = [1244,2051,1301,1465]
           mm(115,1:4) = [2020,2295,3604,2829]
           mm(116,1:4) = [3160,1332,1888,2161]
           mm(117,1:4) = [2785,1832,1836,3365]
           mm(118,1:4) = [2772,2405,1990,361]
           mm(119,1:4) = [1217,3638,2058,2685]
           mm(120,1:4) = [1822,3661,692,3745]
           mm(121,1:4) = [1245,327,1194,2325]
           mm(122,1:4) = [2252,3660,20,3609]
           mm(123,1:4) = [3904,716,3285,3821]
           mm(124,1:4) = [2774,1842,2046,3537]
           mm(125,1:4) = [997,3987,2107,517]
           mm(126,1:4) = [2573,1368,3508,3017]
           mm(127,1:4) = [1148,1848,3525,2141]
           mm(128,1:4) = [545,2366,3801,1537]
           ! Executable Statements
           i1 = iseed(1)
           i2 = iseed(2)
           i3 = iseed(3)
           i4 = iseed(4)
           loop_10: do i = 1,min(n,lv)
           20 continue
              ! multiply the seed by i-th power of the multiplier modulo 2**48
              it4 = i4*mm(i,4)
              it3 = it4/ipw2
              it4 = it4 - ipw2*it3
              it3 = it3 + i3*mm(i,4) + i4*mm(i,3)
              it2 = it3/ipw2
              it3 = it3 - ipw2*it2
              it2 = it2 + i2*mm(i,4) + i3*mm(i,3) + i4*mm(i,2)
              it1 = it2/ipw2
              it2 = it2 - ipw2*it1
              it1 = it1 + i1*mm(i,4) + i2*mm(i,3) + i3*mm(i,2) + i4*mm(i,1)
              it1 = mod(it1,ipw2)
              ! convert 48-bit integer to a realnumber in the interval (0,1,KIND=sp)
              x(i) = r*(real(it1,KIND=sp) + r*(real(it2,KIND=sp) + r*(real(it3,KIND=sp) + &
                        r*real(it4,KIND=sp))))
              if (x(i) == 1.0_sp) then
                 ! if a real number has n bits of precision, and the first
                 ! n bits of the 48-bit integer above happen to be all 1 (which
                 ! will occur about once every 2**n calls), then x( i ) will
                 ! be rounded to exactly one. in ieee single precision arithmetic,
                 ! this will happen relatively often since n = 24.
                 ! since x( i ) is not supposed to return exactly 0.0_sp or 1.0_sp,
                 ! the statistically correct thing to do in this situation is
                 ! simply to iterate again.
                 ! n.b. the case x( i ) = 0.0_sp should not be possible.
                 i1 = i1 + 2
                 i2 = i2 + 2
                 i3 = i3 + 2
                 i4 = i4 + 2
                 goto 20
              end if
           end do loop_10
           ! return final value of seed
           iseed(1) = it1
           iseed(2) = it2
           iseed(3) = it3
           iseed(4) = it4
           return
     end subroutine la_slaruv
     !> DLARUV: returns a vector of n random real numbers from a uniform (0,1)
     !> distribution (n <= 128).
     !> This is an auxiliary routine called by DLARNV and ZLARNV.

     pure subroutine la_dlaruv(iseed,n,x)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(dp),intent(out) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           integer(ilp),parameter :: ipw2 = 4096
           real(dp),parameter :: r = one/ipw2

           ! Local Scalars
           integer(ilp) :: i,i1,i2,i3,i4,it1,it2,it3,it4,j
           ! Local Arrays
           integer(ilp) :: mm(lv,4)
           ! Intrinsic Functions
           intrinsic :: real,min,mod
           ! Data Statements
           mm(1,1:4) = [494,322,2508,2549]
           mm(2,1:4) = [2637,789,3754,1145]
           mm(3,1:4) = [255,1440,1766,2253]
           mm(4,1:4) = [2008,752,3572,305]
           mm(5,1:4) = [1253,2859,2893,3301]
           mm(6,1:4) = [3344,123,307,1065]
           mm(7,1:4) = [4084,1848,1297,3133]
           mm(8,1:4) = [1739,643,3966,2913]
           mm(9,1:4) = [3143,2405,758,3285]
           mm(10,1:4) = [3468,2638,2598,1241]
           mm(11,1:4) = [688,2344,3406,1197]
           mm(12,1:4) = [1657,46,2922,3729]
           mm(13,1:4) = [1238,3814,1038,2501]
           mm(14,1:4) = [3166,913,2934,1673]
           mm(15,1:4) = [1292,3649,2091,541]
           mm(16,1:4) = [3422,339,2451,2753]
           mm(17,1:4) = [1270,3808,1580,949]
           mm(18,1:4) = [2016,822,1958,2361]
           mm(19,1:4) = [154,2832,2055,1165]
           mm(20,1:4) = [2862,3078,1507,4081]
           mm(21,1:4) = [697,3633,1078,2725]
           mm(22,1:4) = [1706,2970,3273,3305]
           mm(23,1:4) = [491,637,17,3069]
           mm(24,1:4) = [931,2249,854,3617]
           mm(25,1:4) = [1444,2081,2916,3733]
           mm(26,1:4) = [444,4019,3971,409]
           mm(27,1:4) = [3577,1478,2889,2157]
           mm(28,1:4) = [3944,242,3831,1361]
           mm(29,1:4) = [2184,481,2621,3973]
           mm(30,1:4) = [1661,2075,1541,1865]
           mm(31,1:4) = [3482,4058,893,2525]
           mm(32,1:4) = [657,622,736,1409]
           mm(33,1:4) = [3023,3376,3992,3445]
           mm(34,1:4) = [3618,812,787,3577]
           mm(35,1:4) = [1267,234,2125,77]
           mm(36,1:4) = [1828,641,2364,3761]
           mm(37,1:4) = [164,4005,2460,2149]
           mm(38,1:4) = [3798,1122,257,1449]
           mm(39,1:4) = [3087,3135,1574,3005]
           mm(40,1:4) = [2400,2640,3912,225]
           mm(41,1:4) = [2870,2302,1216,85]
           mm(42,1:4) = [3876,40,3248,3673]
           mm(43,1:4) = [1905,1832,3401,3117]
           mm(44,1:4) = [1593,2247,2124,3089]
           mm(45,1:4) = [1797,2034,2762,1349]
           mm(46,1:4) = [1234,2637,149,2057]
           mm(47,1:4) = [3460,1287,2245,413]
           mm(48,1:4) = [328,1691,166,65]
           mm(49,1:4) = [2861,496,466,1845]
           mm(50,1:4) = [1950,1597,4018,697]
           mm(51,1:4) = [617,2394,1399,3085]
           mm(52,1:4) = [2070,2584,190,3441]
           mm(53,1:4) = [3331,1843,2879,1573]
           mm(54,1:4) = [769,336,153,3689]
           mm(55,1:4) = [1558,1472,2320,2941]
           mm(56,1:4) = [2412,2407,18,929]
           mm(57,1:4) = [2800,433,712,533]
           mm(58,1:4) = [189,2096,2159,2841]
           mm(59,1:4) = [287,1761,2318,4077]
           mm(60,1:4) = [2045,2810,2091,721]
           mm(61,1:4) = [1227,566,3443,2821]
           mm(62,1:4) = [2838,442,1510,2249]
           mm(63,1:4) = [209,41,449,2397]
           mm(64,1:4) = [2770,1238,1956,2817]
           mm(65,1:4) = [3654,1086,2201,245]
           mm(66,1:4) = [3993,603,3137,1913]
           mm(67,1:4) = [192,840,3399,1997]
           mm(68,1:4) = [2253,3168,1321,3121]
           mm(69,1:4) = [3491,1499,2271,997]
           mm(70,1:4) = [2889,1084,3667,1833]
           mm(71,1:4) = [2857,3438,2703,2877]
           mm(72,1:4) = [2094,2408,629,1633]
           mm(73,1:4) = [1818,1589,2365,981]
           mm(74,1:4) = [688,2391,2431,2009]
           mm(75,1:4) = [1407,288,1113,941]
           mm(76,1:4) = [634,26,3922,2449]
           mm(77,1:4) = [3231,512,2554,197]
           mm(78,1:4) = [815,1456,184,2441]
           mm(79,1:4) = [3524,171,2099,285]
           mm(80,1:4) = [1914,1677,3228,1473]
           mm(81,1:4) = [516,2657,4012,2741]
           mm(82,1:4) = [164,2270,1921,3129]
           mm(83,1:4) = [303,2587,3452,909]
           mm(84,1:4) = [2144,2961,3901,2801]
           mm(85,1:4) = [3480,1970,572,421]
           mm(86,1:4) = [119,1817,3309,4073]
           mm(87,1:4) = [3357,676,3171,2813]
           mm(88,1:4) = [837,1410,817,2337]
           mm(89,1:4) = [2826,3723,3039,1429]
           mm(90,1:4) = [2332,2803,1696,1177]
           mm(91,1:4) = [2089,3185,1256,1901]
           mm(92,1:4) = [3780,184,3715,81]
           mm(93,1:4) = [1700,663,2077,1669]
           mm(94,1:4) = [3712,499,3019,2633]
           mm(95,1:4) = [150,3784,1497,2269]
           mm(96,1:4) = [2000,1631,1101,129]
           mm(97,1:4) = [3375,1925,717,1141]
           mm(98,1:4) = [1621,3912,51,249]
           mm(99,1:4) = [3090,1398,981,3917]
           mm(100,1:4) = [3765,1349,1978,2481]
           mm(101,1:4) = [1149,1441,1813,3941]
           mm(102,1:4) = [3146,2224,3881,2217]
           mm(103,1:4) = [33,2411,76,2749]
           mm(104,1:4) = [3082,1907,3846,3041]
           mm(105,1:4) = [2741,3192,3694,1877]
           mm(106,1:4) = [359,2786,1682,345]
           mm(107,1:4) = [3316,382,124,2861]
           mm(108,1:4) = [1749,37,1660,1809]
           mm(109,1:4) = [185,759,3997,3141]
           mm(110,1:4) = [2784,2948,479,2825]
           mm(111,1:4) = [2202,1862,1141,157]
           mm(112,1:4) = [2199,3802,886,2881]
           mm(113,1:4) = [1364,2423,3514,3637]
           mm(114,1:4) = [1244,2051,1301,1465]
           mm(115,1:4) = [2020,2295,3604,2829]
           mm(116,1:4) = [3160,1332,1888,2161]
           mm(117,1:4) = [2785,1832,1836,3365]
           mm(118,1:4) = [2772,2405,1990,361]
           mm(119,1:4) = [1217,3638,2058,2685]
           mm(120,1:4) = [1822,3661,692,3745]
           mm(121,1:4) = [1245,327,1194,2325]
           mm(122,1:4) = [2252,3660,20,3609]
           mm(123,1:4) = [3904,716,3285,3821]
           mm(124,1:4) = [2774,1842,2046,3537]
           mm(125,1:4) = [997,3987,2107,517]
           mm(126,1:4) = [2573,1368,3508,3017]
           mm(127,1:4) = [1148,1848,3525,2141]
           mm(128,1:4) = [545,2366,3801,1537]
           ! Executable Statements
           i1 = iseed(1)
           i2 = iseed(2)
           i3 = iseed(3)
           i4 = iseed(4)
           loop_10: do i = 1,min(n,lv)
           20 continue
              ! multiply the seed by i-th power of the multiplier modulo 2**48
              it4 = i4*mm(i,4)
              it3 = it4/ipw2
              it4 = it4 - ipw2*it3
              it3 = it3 + i3*mm(i,4) + i4*mm(i,3)
              it2 = it3/ipw2
              it3 = it3 - ipw2*it2
              it2 = it2 + i2*mm(i,4) + i3*mm(i,3) + i4*mm(i,2)
              it1 = it2/ipw2
              it2 = it2 - ipw2*it1
              it1 = it1 + i1*mm(i,4) + i2*mm(i,3) + i3*mm(i,2) + i4*mm(i,1)
              it1 = mod(it1,ipw2)
              ! convert 48-bit integer to a realnumber in the interval (0,1,KIND=dp)
              x(i) = r*(real(it1,KIND=dp) + r*(real(it2,KIND=dp) + r*(real(it3,KIND=dp) + &
                        r*real(it4,KIND=dp))))
              if (x(i) == 1.0_dp) then
                 ! if a real number has n bits of precision, and the first
                 ! n bits of the 48-bit integer above happen to be all 1 (which
                 ! will occur about once every 2**n calls), then x( i ) will
                 ! be rounded to exactly one.
                 ! since x( i ) is not supposed to return exactly 0.0_dp or 1.0_dp,
                 ! the statistically correct thing to do in this situation is
                 ! simply to iterate again.
                 ! n.b. the case x( i ) = 0.0_dp should not be possible.
                 i1 = i1 + 2
                 i2 = i2 + 2
                 i3 = i3 + 2
                 i4 = i4 + 2
                 goto 20
              end if
           end do loop_10
           ! return final value of seed
           iseed(1) = it1
           iseed(2) = it2
           iseed(3) = it3
           iseed(4) = it4
           return
     end subroutine la_dlaruv
     !> QLARUV: returns a vector of n random real numbers from a uniform (0,1)
     !> distribution (n <= 128).
     !> This is an auxiliary routine called by QLARNV and WLARNV.

     pure subroutine la_qlaruv(iseed,n,x)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(qp),intent(out) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           integer(ilp),parameter :: ipw2 = 4096
           real(qp),parameter :: r = one/ipw2

           ! Local Scalars
           integer(ilp) :: i,i1,i2,i3,i4,it1,it2,it3,it4,j
           ! Local Arrays
           integer(ilp) :: mm(lv,4)
           ! Intrinsic Functions
           intrinsic :: real,min,mod
           ! Data Statements
           mm(1,1:4) = [494,322,2508,2549]
           mm(2,1:4) = [2637,789,3754,1145]
           mm(3,1:4) = [255,1440,1766,2253]
           mm(4,1:4) = [2008,752,3572,305]
           mm(5,1:4) = [1253,2859,2893,3301]
           mm(6,1:4) = [3344,123,307,1065]
           mm(7,1:4) = [4084,1848,1297,3133]
           mm(8,1:4) = [1739,643,3966,2913]
           mm(9,1:4) = [3143,2405,758,3285]
           mm(10,1:4) = [3468,2638,2598,1241]
           mm(11,1:4) = [688,2344,3406,1197]
           mm(12,1:4) = [1657,46,2922,3729]
           mm(13,1:4) = [1238,3814,1038,2501]
           mm(14,1:4) = [3166,913,2934,1673]
           mm(15,1:4) = [1292,3649,2091,541]
           mm(16,1:4) = [3422,339,2451,2753]
           mm(17,1:4) = [1270,3808,1580,949]
           mm(18,1:4) = [2016,822,1958,2361]
           mm(19,1:4) = [154,2832,2055,1165]
           mm(20,1:4) = [2862,3078,1507,4081]
           mm(21,1:4) = [697,3633,1078,2725]
           mm(22,1:4) = [1706,2970,3273,3305]
           mm(23,1:4) = [491,637,17,3069]
           mm(24,1:4) = [931,2249,854,3617]
           mm(25,1:4) = [1444,2081,2916,3733]
           mm(26,1:4) = [444,4019,3971,409]
           mm(27,1:4) = [3577,1478,2889,2157]
           mm(28,1:4) = [3944,242,3831,1361]
           mm(29,1:4) = [2184,481,2621,3973]
           mm(30,1:4) = [1661,2075,1541,1865]
           mm(31,1:4) = [3482,4058,893,2525]
           mm(32,1:4) = [657,622,736,1409]
           mm(33,1:4) = [3023,3376,3992,3445]
           mm(34,1:4) = [3618,812,787,3577]
           mm(35,1:4) = [1267,234,2125,77]
           mm(36,1:4) = [1828,641,2364,3761]
           mm(37,1:4) = [164,4005,2460,2149]
           mm(38,1:4) = [3798,1122,257,1449]
           mm(39,1:4) = [3087,3135,1574,3005]
           mm(40,1:4) = [2400,2640,3912,225]
           mm(41,1:4) = [2870,2302,1216,85]
           mm(42,1:4) = [3876,40,3248,3673]
           mm(43,1:4) = [1905,1832,3401,3117]
           mm(44,1:4) = [1593,2247,2124,3089]
           mm(45,1:4) = [1797,2034,2762,1349]
           mm(46,1:4) = [1234,2637,149,2057]
           mm(47,1:4) = [3460,1287,2245,413]
           mm(48,1:4) = [328,1691,166,65]
           mm(49,1:4) = [2861,496,466,1845]
           mm(50,1:4) = [1950,1597,4018,697]
           mm(51,1:4) = [617,2394,1399,3085]
           mm(52,1:4) = [2070,2584,190,3441]
           mm(53,1:4) = [3331,1843,2879,1573]
           mm(54,1:4) = [769,336,153,3689]
           mm(55,1:4) = [1558,1472,2320,2941]
           mm(56,1:4) = [2412,2407,18,929]
           mm(57,1:4) = [2800,433,712,533]
           mm(58,1:4) = [189,2096,2159,2841]
           mm(59,1:4) = [287,1761,2318,4077]
           mm(60,1:4) = [2045,2810,2091,721]
           mm(61,1:4) = [1227,566,3443,2821]
           mm(62,1:4) = [2838,442,1510,2249]
           mm(63,1:4) = [209,41,449,2397]
           mm(64,1:4) = [2770,1238,1956,2817]
           mm(65,1:4) = [3654,1086,2201,245]
           mm(66,1:4) = [3993,603,3137,1913]
           mm(67,1:4) = [192,840,3399,1997]
           mm(68,1:4) = [2253,3168,1321,3121]
           mm(69,1:4) = [3491,1499,2271,997]
           mm(70,1:4) = [2889,1084,3667,1833]
           mm(71,1:4) = [2857,3438,2703,2877]
           mm(72,1:4) = [2094,2408,629,1633]
           mm(73,1:4) = [1818,1589,2365,981]
           mm(74,1:4) = [688,2391,2431,2009]
           mm(75,1:4) = [1407,288,1113,941]
           mm(76,1:4) = [634,26,3922,2449]
           mm(77,1:4) = [3231,512,2554,197]
           mm(78,1:4) = [815,1456,184,2441]
           mm(79,1:4) = [3524,171,2099,285]
           mm(80,1:4) = [1914,1677,3228,1473]
           mm(81,1:4) = [516,2657,4012,2741]
           mm(82,1:4) = [164,2270,1921,3129]
           mm(83,1:4) = [303,2587,3452,909]
           mm(84,1:4) = [2144,2961,3901,2801]
           mm(85,1:4) = [3480,1970,572,421]
           mm(86,1:4) = [119,1817,3309,4073]
           mm(87,1:4) = [3357,676,3171,2813]
           mm(88,1:4) = [837,1410,817,2337]
           mm(89,1:4) = [2826,3723,3039,1429]
           mm(90,1:4) = [2332,2803,1696,1177]
           mm(91,1:4) = [2089,3185,1256,1901]
           mm(92,1:4) = [3780,184,3715,81]
           mm(93,1:4) = [1700,663,2077,1669]
           mm(94,1:4) = [3712,499,3019,2633]
           mm(95,1:4) = [150,3784,1497,2269]
           mm(96,1:4) = [2000,1631,1101,129]
           mm(97,1:4) = [3375,1925,717,1141]
           mm(98,1:4) = [1621,3912,51,249]
           mm(99,1:4) = [3090,1398,981,3917]
           mm(100,1:4) = [3765,1349,1978,2481]
           mm(101,1:4) = [1149,1441,1813,3941]
           mm(102,1:4) = [3146,2224,3881,2217]
           mm(103,1:4) = [33,2411,76,2749]
           mm(104,1:4) = [3082,1907,3846,3041]
           mm(105,1:4) = [2741,3192,3694,1877]
           mm(106,1:4) = [359,2786,1682,345]
           mm(107,1:4) = [3316,382,124,2861]
           mm(108,1:4) = [1749,37,1660,1809]
           mm(109,1:4) = [185,759,3997,3141]
           mm(110,1:4) = [2784,2948,479,2825]
           mm(111,1:4) = [2202,1862,1141,157]
           mm(112,1:4) = [2199,3802,886,2881]
           mm(113,1:4) = [1364,2423,3514,3637]
           mm(114,1:4) = [1244,2051,1301,1465]
           mm(115,1:4) = [2020,2295,3604,2829]
           mm(116,1:4) = [3160,1332,1888,2161]
           mm(117,1:4) = [2785,1832,1836,3365]
           mm(118,1:4) = [2772,2405,1990,361]
           mm(119,1:4) = [1217,3638,2058,2685]
           mm(120,1:4) = [1822,3661,692,3745]
           mm(121,1:4) = [1245,327,1194,2325]
           mm(122,1:4) = [2252,3660,20,3609]
           mm(123,1:4) = [3904,716,3285,3821]
           mm(124,1:4) = [2774,1842,2046,3537]
           mm(125,1:4) = [997,3987,2107,517]
           mm(126,1:4) = [2573,1368,3508,3017]
           mm(127,1:4) = [1148,1848,3525,2141]
           mm(128,1:4) = [545,2366,3801,1537]
           ! Executable Statements
           i1 = iseed(1)
           i2 = iseed(2)
           i3 = iseed(3)
           i4 = iseed(4)
           loop_10: do i = 1,min(n,lv)
           20 continue
              ! multiply the seed by i-th power of the multiplier modulo 2**48
              it4 = i4*mm(i,4)
              it3 = it4/ipw2
              it4 = it4 - ipw2*it3
              it3 = it3 + i3*mm(i,4) + i4*mm(i,3)
              it2 = it3/ipw2
              it3 = it3 - ipw2*it2
              it2 = it2 + i2*mm(i,4) + i3*mm(i,3) + i4*mm(i,2)
              it1 = it2/ipw2
              it2 = it2 - ipw2*it1
              it1 = it1 + i1*mm(i,4) + i2*mm(i,3) + i3*mm(i,2) + i4*mm(i,1)
              it1 = mod(it1,ipw2)
              ! convert 48-bit integer to a realnumber in the interval (0,1,KIND=qp)
              x(i) = r*(real(it1,KIND=qp) + r*(real(it2,KIND=qp) + r*(real(it3,KIND=qp) + &
                        r*real(it4,KIND=qp))))
              if (x(i) == 1.0_qp) then
                 ! if a real number has n bits of precision, and the first
                 ! n bits of the 48-bit integer above happen to be all 1 (which
                 ! will occur about once every 2**n calls), then x( i ) will
                 ! be rounded to exactly one.
                 ! since x( i ) is not supposed to return exactly 0.0_qp or 1.0_qp,
                 ! the statistically correct thing to do in this situation is
                 ! simply to iterate again.
                 ! n.b. the case x( i ) = 0.0_qp should not be possible.
                 i1 = i1 + 2
                 i2 = i2 + 2
                 i3 = i3 + 2
                 i4 = i4 + 2
                 goto 20
              end if
           end do loop_10
           ! return final value of seed
           iseed(1) = it1
           iseed(2) = it2
           iseed(3) = it3
           iseed(4) = it4
           return
     end subroutine la_qlaruv

     !> SLASET: initializes an m-by-n matrix A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_slaset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           real(sp),intent(in) :: alpha,beta
           ! Array Arguments
           real(sp),intent(out) :: a(lda,*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the strictly upper triangular or trapezoidal part of the
              ! array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the strictly lower triangular or trapezoidal part of the
              ! array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
           else
              ! set the leading m-by-n submatrix to alpha.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
           end if
           ! set the first min(m,n) diagonal elements to beta.
           do i = 1,min(m,n)
              a(i,i) = beta
           end do
           return
     end subroutine la_slaset
     !> DLASET: initializes an m-by-n matrix A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_dlaset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           real(dp),intent(in) :: alpha,beta
           ! Array Arguments
           real(dp),intent(out) :: a(lda,*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the strictly upper triangular or trapezoidal part of the
              ! array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the strictly lower triangular or trapezoidal part of the
              ! array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
           else
              ! set the leading m-by-n submatrix to alpha.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
           end if
           ! set the first min(m,n) diagonal elements to beta.
           do i = 1,min(m,n)
              a(i,i) = beta
           end do
           return
     end subroutine la_dlaset
     !> QLASET: initializes an m-by-n matrix A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_qlaset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           real(qp),intent(in) :: alpha,beta
           ! Array Arguments
           real(qp),intent(out) :: a(lda,*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the strictly upper triangular or trapezoidal part of the
              ! array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the strictly lower triangular or trapezoidal part of the
              ! array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
           else
              ! set the leading m-by-n submatrix to alpha.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
           end if
           ! set the first min(m,n) diagonal elements to beta.
           do i = 1,min(m,n)
              a(i,i) = beta
           end do
           return
     end subroutine la_qlaset

     !> DLAT2S: converts a DOUBLE PRECISION triangular matrix, SA, to a SINGLE
     !> PRECISION triangular matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> DLAS2S checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_dlat2s(uplo,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,n
           ! Array Arguments
           real(sp),intent(out) :: sa(ldsa,*)
           real(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: rmax
           logical(lk) :: upper
           ! Executable Statements
           rmax = la_slamch('O')
           upper = la_lsame(uplo,'U')
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           end if
           50 continue
           return
     end subroutine la_dlat2s
     !> QLAT2D: converts a QUAD PRECISION triangular matrix, SA, to a SINGLE
     !> PRECISION triangular matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> DLAS2S checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_qlat2d(uplo,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,n
           ! Array Arguments
           real(dp),intent(out) :: sa(ldsa,*)
           real(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: rmax
           logical(lk) :: upper
           ! Executable Statements
           rmax = la_dlamch('O')
           upper = la_lsame(uplo,'U')
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    if ((a(i,j) < -rmax) .or. (a(i,j) > rmax)) then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           end if
           50 continue
           return
     end subroutine la_qlat2d

     !> STFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_stfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(out) :: ap(0:*)
           real(sp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('STFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_stfttp
     !> DTFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_dtfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(out) :: ap(0:*)
           real(dp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('DTFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dtfttp
     !> QTFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_qtfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(out) :: ap(0:*)
           real(qp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('QTFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qtfttp

     !> STFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_stfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(sp),intent(out) :: a(0:lda - 1,0:*)
           real(sp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('STFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 a(0,0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_stfttr
     !> DTFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_dtfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(dp),intent(out) :: a(0:lda - 1,0:*)
           real(dp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DTFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 a(0,0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dtfttr
     !> QTFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_qtfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(qp),intent(out) :: a(0:lda - 1,0:*)
           real(qp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QTFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 a(0,0) = arf(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = arf(ij)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qtfttr

     !> STPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_stpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(in) :: ap(0:*)
           real(sp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('STPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = ap(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_stpttf
     !> DTPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_dtpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(in) :: ap(0:*)
           real(dp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('DTPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = ap(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dtpttf
     !> QTPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_qtpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(in) :: ap(0:*)
           real(qp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('QTPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = ap(0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qtpttf

     !> STPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_stpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(sp),intent(out) :: a(lda,*)
           real(sp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('STPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_stpttr
     !> DTPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_dtpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(dp),intent(out) :: a(lda,*)
           real(dp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DTPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_dtpttr
     !> QTPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_qtpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(qp),intent(out) :: a(lda,*)
           real(qp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QTPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_qtpttr

     !> STRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_strttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(sp),intent(in) :: a(0:lda - 1,0:*)
           real(sp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('STRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 arf(0) = a(0,0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = a(n2 + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = a(j - n1,l)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = a(n2 + j,l)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = a(k + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = a(j - k,l)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = a(k + 1 + j,l)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_strttf
     !> DTRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_dtrttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(dp),intent(in) :: a(0:lda - 1,0:*)
           real(dp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DTRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 arf(0) = a(0,0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = a(n2 + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = a(j - n1,l)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = a(n2 + j,l)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = a(k + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = a(j - k,l)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = a(k + 1 + j,l)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dtrttf
     !> QTRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_qtrttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(qp),intent(in) :: a(0:lda - 1,0:*)
           real(qp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QTRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 arf(0) = a(0,0)
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = a(n2 + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = a(j - n1,l)
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = a(n2 + j,l)
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = a(k + j,i)
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = a(j - k,l)
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = a(j,i)
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = a(k + 1 + j,l)
                          ij = ij + 1
                       end do
                    end do
                    ! note that here, on exit of the loop, j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qtrttf

     !> STRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_strttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('STRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_strttp
     !> DTRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_dtrttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DTRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_dtrttp
     !> QTRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_qtrttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QTRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_qtrttp

     !> SLAG2D: converts a SINGLE PRECISION matrix, SA, to a DOUBLE
     !> PRECISION matrix, A.
     !> Note that while it is possible to overflow while converting
     !> from double to single, it is not possible to overflow when
     !> converting from single to double.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_slag2d(m,n,sa,ldsa,a,lda,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           real(sp),intent(in) :: sa(ldsa,*)
           real(dp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           info = 0
           do j = 1,n
              do i = 1,m
                 a(i,j) = sa(i,j)
              end do
           end do
           return
     end subroutine la_slag2d
     !> DLAG2Q: converts a SINGLE PRECISION matrix, SA, to a DOUBLE
     !> PRECISION matrix, A.
     !> Note that while it is possible to overflow while converting
     !> from double to single, it is not possible to overflow when
     !> converting from single to double.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_dlag2q(m,n,sa,ldsa,a,lda,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           real(dp),intent(in) :: sa(ldsa,*)
           real(qp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           info = 0
           do j = 1,n
              do i = 1,m
                 a(i,j) = sa(i,j)
              end do
           end do
           return
     end subroutine la_dlag2q

     !> SLARNV: returns a vector of n random real numbers from a uniform or
     !> normal distribution.

     pure subroutine la_slarnv(idist,iseed,n,x)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(sp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(sp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_sp

           ! Local Scalars
           integer(ilp) :: i,il,il2,iv
           ! Local Arrays
           real(sp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cos,log,min,sqrt
           ! Executable Statements
           do 40 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              if (idist == 3) then
                 il2 = 2*il
              else
                 il2 = il
              end if
              ! call la_slaruv to generate il2 numbers from a uniform (0,1)
              ! distribution (il2 <= lv)
              call la_slaruv(iseed,il2,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = u(i)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = two*u(i) - one
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*cos(twopi*u(2*i))
                 end do
              end if
              40 continue
           return
     end subroutine la_slarnv
     !> DLARNV: returns a vector of n random real numbers from a uniform or
     !> normal distribution.

     pure subroutine la_dlarnv(idist,iseed,n,x)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(dp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(dp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_dp

           ! Local Scalars
           integer(ilp) :: i,il,il2,iv
           ! Local Arrays
           real(dp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cos,log,min,sqrt
           ! Executable Statements
           do 40 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              if (idist == 3) then
                 il2 = 2*il
              else
                 il2 = il
              end if
              ! call la_dlaruv to generate il2 numbers from a uniform (0,1)
              ! distribution (il2 <= lv)
              call la_dlaruv(iseed,il2,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = u(i)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = two*u(i) - one
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*cos(twopi*u(2*i))
                 end do
              end if
              40 continue
           return
     end subroutine la_dlarnv
     !> QLARNV: returns a vector of n random real numbers from a uniform or
     !> normal distribution.

     pure subroutine la_qlarnv(idist,iseed,n,x)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           real(qp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(qp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_qp

           ! Local Scalars
           integer(ilp) :: i,il,il2,iv
           ! Local Arrays
           real(qp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cos,log,min,sqrt
           ! Executable Statements
           do 40 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              if (idist == 3) then
                 il2 = 2*il
              else
                 il2 = il
              end if
              ! call la_qlaruv to generate il2 numbers from a uniform (0,1)
              ! distribution (il2 <= lv)
              call la_qlaruv(iseed,il2,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = u(i)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = two*u(i) - one
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*cos(twopi*u(2*i))
                 end do
              end if
              40 continue
           return
     end subroutine la_qlarnv

     !> CLAG2Z: converts a COMPLEX matrix, SA, to a COMPLEX*16 matrix, A.
     !> Note that while it is possible to overflow while converting
     !> from double to single, it is not possible to overflow when
     !> converting from single to double.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_clag2z(m,n,sa,ldsa,a,lda,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           complex(sp),intent(in) :: sa(ldsa,*)
           complex(dp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           info = 0
           do j = 1,n
              do i = 1,m
                 a(i,j) = sa(i,j)
              end do
           end do
           return
     end subroutine la_clag2z
     !> ZLAG2W: converts a COMPLEX matrix, SA, to a COMPLEX*16 matrix, A.
     !> Note that while it is possible to overflow while converting
     !> from double to single, it is not possible to overflow when
     !> converting from single to double.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_zlag2w(m,n,sa,ldsa,a,lda,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           complex(dp),intent(in) :: sa(ldsa,*)
           complex(qp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           info = 0
           do j = 1,n
              do i = 1,m
                 a(i,j) = sa(i,j)
              end do
           end do
           return
     end subroutine la_zlag2w

     !> CLACP2: copies all or part of a real two-dimensional matrix A to a
     !> complex matrix B.

     pure subroutine la_clacp2(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_clacp2
     !> ZLACP2: copies all or part of a real two-dimensional matrix A to a
     !> complex matrix B.

     pure subroutine la_zlacp2(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_zlacp2
     !> WLACP2: copies all or part of a real two-dimensional matrix A to a
     !> complex matrix B.

     pure subroutine la_wlacp2(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_wlacp2

     !> CLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_clacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_clacpy
     !> ZLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_zlacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_zlacpy
     !> WLACPY: copies all or part of a two-dimensional matrix A to another
     !> matrix B.

     pure subroutine la_wlacpy(uplo,m,n,a,lda,b,ldb)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,ldb,m,n
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: b(ldb,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              do j = 1,n
                 do i = 1,min(j,m)
                    b(i,j) = a(i,j)
                 end do
              end do
           else if (la_lsame(uplo,'L')) then
              do j = 1,n
                 do i = j,m
                    b(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = 1,m
                    b(i,j) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_wlacpy

     !> ZLAG2C: converts a COMPLEX*16 matrix, SA, to a COMPLEX matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> ZLAG2C checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_zlag2c(m,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           complex(sp),intent(out) :: sa(ldsa,*)
           complex(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: rmax
           ! Intrinsic Functions
           intrinsic :: real,aimag
           ! Executable Statements
           rmax = la_slamch('O')
           do j = 1,n
              do i = 1,m
                 if ((real(a(i,j),KIND=dp) < -rmax) .or. (real(a(i,j),KIND=dp) > rmax) &
                           .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) then
                    info = 1
                    go to 30
                 end if
                 sa(i,j) = a(i,j)
              end do
           end do
           info = 0
           30 continue
           return
     end subroutine la_zlag2c
     !> WLAG2Z: converts a COMPLEX*16 matrix, SA, to a COMPLEX matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> WLAG2Z checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_wlag2z(m,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,m,n
           ! Array Arguments
           complex(dp),intent(out) :: sa(ldsa,*)
           complex(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: rmax
           ! Intrinsic Functions
           intrinsic :: real,aimag
           ! Executable Statements
           rmax = la_dlamch('O')
           do j = 1,n
              do i = 1,m
                 if ((real(a(i,j),KIND=qp) < -rmax) .or. (real(a(i,j),KIND=qp) > rmax) &
                           .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) then
                    info = 1
                    go to 30
                 end if
                 sa(i,j) = a(i,j)
              end do
           end do
           info = 0
           30 continue
           return
     end subroutine la_wlag2z

     !> CLARNV: returns a vector of n random complex numbers from a uniform or
     !> normal distribution.

     pure subroutine la_clarnv(idist,iseed,n,x)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           complex(sp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(sp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_sp

           ! Local Scalars
           integer(ilp) :: i,il,iv
           ! Local Arrays
           real(sp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cmplx,exp,log,min,sqrt
           ! Executable Statements
           do 60 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              ! call la_slaruv to generate 2*il realnumbers from a uniform (0,1,KIND=sp)
              ! distribution (2*il <= lv)
              call la_slaruv(iseed,2*il,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = cmplx(u(2*i - 1),u(2*i),KIND=sp)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = cmplx(two*u(2*i - 1) - one,two*u(2*i) - one,KIND=sp)
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*exp(cmplx(zero,twopi*u(2*i), &
                              KIND=sp))
                 end do
              else if (idist == 4) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit disk
                 do i = 1,il
                    x(iv + i - 1) = sqrt(u(2*i - 1))*exp(cmplx(zero,twopi*u(2*i),KIND=sp))

                 end do
              else if (idist == 5) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit circle
                 do i = 1,il
                    x(iv + i - 1) = exp(cmplx(zero,twopi*u(2*i),KIND=sp))
                 end do
              end if
              60 continue
           return
     end subroutine la_clarnv
     !> ZLARNV: returns a vector of n random complex numbers from a uniform or
     !> normal distribution.

     pure subroutine la_zlarnv(idist,iseed,n,x)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           complex(dp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(dp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_dp

           ! Local Scalars
           integer(ilp) :: i,il,iv
           ! Local Arrays
           real(dp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cmplx,exp,log,min,sqrt
           ! Executable Statements
           do 60 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              ! call la_dlaruv to generate 2*il realnumbers from a uniform (0,1,KIND=dp)
              ! distribution (2*il <= lv)
              call la_dlaruv(iseed,2*il,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = cmplx(u(2*i - 1),u(2*i),KIND=dp)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = cmplx(two*u(2*i - 1) - one,two*u(2*i) - one,KIND=dp)
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*exp(cmplx(zero,twopi*u(2*i), &
                              KIND=dp))
                 end do
              else if (idist == 4) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit disk
                 do i = 1,il
                    x(iv + i - 1) = sqrt(u(2*i - 1))*exp(cmplx(zero,twopi*u(2*i),KIND=dp))

                 end do
              else if (idist == 5) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit circle
                 do i = 1,il
                    x(iv + i - 1) = exp(cmplx(zero,twopi*u(2*i),KIND=dp))
                 end do
              end if
              60 continue
           return
     end subroutine la_zlarnv
     !> WLARNV: returns a vector of n random complex numbers from a uniform or
     !> normal distribution.

     pure subroutine la_wlarnv(idist,iseed,n,x)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: idist,n
           ! Array Arguments
           integer(ilp),intent(inout) :: iseed(4)
           complex(qp),intent(out) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: lv = 128
           real(qp),parameter :: twopi = 6.28318530717958647692528676655900576839e+0_qp

           ! Local Scalars
           integer(ilp) :: i,il,iv
           ! Local Arrays
           real(qp) :: u(lv)
           ! Intrinsic Functions
           intrinsic :: cmplx,exp,log,min,sqrt
           ! Executable Statements
           do 60 iv = 1,n,lv/2
              il = min(lv/2,n - iv + 1)
              ! call la_qlaruv to generate 2*il realnumbers from a uniform (0,1,KIND=qp)
              ! distribution (2*il <= lv)
              call la_qlaruv(iseed,2*il,u)
              if (idist == 1) then
                 ! copy generated numbers
                 do i = 1,il
                    x(iv + i - 1) = cmplx(u(2*i - 1),u(2*i),KIND=qp)
                 end do
              else if (idist == 2) then
                 ! convert generated numbers to uniform (-1,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = cmplx(two*u(2*i - 1) - one,two*u(2*i) - one,KIND=qp)
                 end do
              else if (idist == 3) then
                 ! convert generated numbers to normal (0,1) distribution
                 do i = 1,il
                    x(iv + i - 1) = sqrt(-two*log(u(2*i - 1)))*exp(cmplx(zero,twopi*u(2*i), &
                              KIND=qp))
                 end do
              else if (idist == 4) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit disk
                 do i = 1,il
                    x(iv + i - 1) = sqrt(u(2*i - 1))*exp(cmplx(zero,twopi*u(2*i),KIND=qp))

                 end do
              else if (idist == 5) then
                 ! convert generated numbers to complex numbers uniformly
                 ! distributed on the unit circle
                 do i = 1,il
                    x(iv + i - 1) = exp(cmplx(zero,twopi*u(2*i),KIND=qp))
                 end do
              end if
              60 continue
           return
     end subroutine la_wlarnv

     !> CLASET: initializes a 2-D array A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_claset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           complex(sp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(sp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the diagonal to beta and the strictly upper triangular
              ! part of the array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the diagonal to beta and the strictly lower triangular
              ! part of the array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else
              ! set the array to beta on the diagonal and alpha on the
              ! offdiagonal.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(m,n)
                 a(i,i) = beta
              end do
           end if
           return
     end subroutine la_claset
     !> ZLASET: initializes a 2-D array A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_zlaset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           complex(dp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(dp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the diagonal to beta and the strictly upper triangular
              ! part of the array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the diagonal to beta and the strictly lower triangular
              ! part of the array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else
              ! set the array to beta on the diagonal and alpha on the
              ! offdiagonal.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(m,n)
                 a(i,i) = beta
              end do
           end if
           return
     end subroutine la_zlaset
     !> WLASET: initializes a 2-D array A to BETA on the diagonal and
     !> ALPHA on the offdiagonals.

     pure subroutine la_wlaset(uplo,m,n,alpha,beta,a,lda)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,m,n
           complex(qp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(qp),intent(out) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if (la_lsame(uplo,'U')) then
              ! set the diagonal to beta and the strictly upper triangular
              ! part of the array to alpha.
              do j = 2,n
                 do i = 1,min(j - 1,m)
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else if (la_lsame(uplo,'L')) then
              ! set the diagonal to beta and the strictly lower triangular
              ! part of the array to alpha.
              do j = 1,min(m,n)
                 do i = j + 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(n,m)
                 a(i,i) = beta
              end do
           else
              ! set the array to beta on the diagonal and alpha on the
              ! offdiagonal.
              do j = 1,n
                 do i = 1,m
                    a(i,j) = alpha
                 end do
              end do
              do i = 1,min(m,n)
                 a(i,i) = beta
              end do
           end if
           return
     end subroutine la_wlaset

     !> ZLAT2C: converts a COMPLEX*16 triangular matrix, SA, to a COMPLEX
     !> triangular matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> ZLAT2C checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_zlat2c(uplo,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,n
           ! Array Arguments
           complex(sp),intent(out) :: sa(ldsa,*)
           complex(dp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: rmax
           logical(lk) :: upper
           ! Intrinsic Functions
           intrinsic :: real,aimag
           ! Executable Statements
           rmax = la_slamch('O')
           upper = la_lsame(uplo,'U')
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    if ((real(a(i,j),KIND=dp) < -rmax) .or. (real(a(i,j),KIND=dp) > rmax) &
                    .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) &
                              then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    if ((real(a(i,j),KIND=dp) < -rmax) .or. (real(a(i,j),KIND=dp) > rmax) &
                    .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) &
                              then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           end if
           50 continue
           return
     end subroutine la_zlat2c
     !> WLAT2Z: converts a COMPLEX*16 triangular matrix, SA, to a COMPLEX
     !> triangular matrix, A.
     !> RMAX is the overflow for the SINGLE PRECISION arithmetic
     !> WLAT2Z checks that all the entries of A are between -RMAX and
     !> RMAX. If not the conversion is aborted and a flag is raised.
     !> This is an auxiliary routine so there is no argument checking.

     pure subroutine la_wlat2z(uplo,n,a,lda,sa,ldsa,info)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldsa,n
           ! Array Arguments
           complex(dp),intent(out) :: sa(ldsa,*)
           complex(qp),intent(in) :: a(lda,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: rmax
           logical(lk) :: upper
           ! Intrinsic Functions
           intrinsic :: real,aimag
           ! Executable Statements
           rmax = la_dlamch('O')
           upper = la_lsame(uplo,'U')
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    if ((real(a(i,j),KIND=qp) < -rmax) .or. (real(a(i,j),KIND=qp) > rmax) &
                    .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) &
                              then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    if ((real(a(i,j),KIND=qp) < -rmax) .or. (real(a(i,j),KIND=qp) > rmax) &
                    .or. (aimag(a(i,j)) < -rmax) .or. (aimag(a(i,j)) > rmax)) &
                              then
                       info = 1
                       go to 50
                    end if
                    sa(i,j) = a(i,j)
                 end do
              end do
           end if
           50 continue
           return
     end subroutine la_wlat2z

     !> CTFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_ctfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(out) :: ap(0:*)
           complex(sp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Intrinsic Functions
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('CTFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = conjg(arf(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ctfttp
     !> ZTFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_ztfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(out) :: ap(0:*)
           complex(dp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Intrinsic Functions
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('ZTFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = conjg(arf(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ztfttp
     !> WTFTTP: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard packed format (TP).

     pure subroutine la_wtfttp(transr,uplo,n,arf,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(out) :: ap(0:*)
           complex(qp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Intrinsic Functions
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('WTFTTP',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 ap(0) = arf(0)
              else
                 ap(0) = conjg(arf(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          ap(ijp) = arf(ij)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          ap(ijp) = conjg(arf(ij))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wtfttp

     !> CTFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_ctfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(sp),intent(out) :: a(0:lda - 1,0:*)
           complex(sp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CTFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    a(0,0) = arf(0)
                 else
                    a(0,0) = conjg(arf(0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ctfttr
     !> ZTFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_ztfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(dp),intent(out) :: a(0:lda - 1,0:*)
           complex(dp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZTFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    a(0,0) = arf(0)
                 else
                    a(0,0) = conjg(arf(0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ztfttr
     !> WTFTTR: copies a triangular matrix A from rectangular full packed
     !> format (TF) to standard full format (TR).

     pure subroutine la_wtfttr(transr,uplo,n,arf,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(qp),intent(out) :: a(0:lda - 1,0:*)
           complex(qp),intent(in) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt,nx2,np1x2
           integer(ilp) :: i,j,l,ij
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WTFTTR',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    a(0,0) = arf(0)
                 else
                    a(0,0) = conjg(arf(0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          a(n2 + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          a(j - n1,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          a(i,n1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          a(n2 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          a(k + j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          a(j - k,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          a(i,k + 1 + j) = arf(ij)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          a(j,i) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          a(i,j) = arf(ij)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          a(k + 1 + j,l) = conjg(arf(ij))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       a(i,j) = arf(ij)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wtfttr

     !> CTPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_ctpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(in) :: ap(0:*)
           complex(sp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('CTPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = conjg(ap(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ctpttf
     !> ZTPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_ztpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(in) :: ap(0:*)
           complex(dp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('ZTPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = conjg(ap(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ztpttf
     !> WTPTTF: copies a triangular matrix A from standard packed format (TP)
     !> to rectangular full packed format (TF).

     pure subroutine la_wtpttf(transr,uplo,n,ap,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(in) :: ap(0:*)
           complex(qp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k,nt
           integer(ilp) :: i,j,ij
           integer(ilp) :: ijp,jp,lda,js
           ! Intrinsic Functions
           intrinsic :: conjg,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('WTPTTF',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (normaltransr) then
                 arf(0) = ap(0)
              else
                 arf(0) = conjg(ap(0))
              end if
              return
           end if
           ! size of array arf(0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           ! set lda of arf^c; arf^c is (0:(n+1)/2-1,0:n-noe)
           ! where noe = 0 if n is even, noe = 1 if n is odd
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              lda = n + 1
           else
              nisodd = .true.
              lda = n
           end if
           ! arf^c has lda rows and n+1-noe cols
           if (.not. normaltransr) lda = (n + 1)/2
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda = n
                    ijp = 0
                    jp = 0
                    do j = 0,n2
                       do i = j,n - 1
                          ij = i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,n2 - 1
                       do j = 1 + i,n2
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    ijp = 0
                    do j = 0,n1 - 1
                       ij = n2 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = n1,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ijp = 0
                    do i = 0,n2
                       do ij = i*(lda + 1),n*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 1
                    do j = 0,n2 - 1
                       do ij = js,js + n2 - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda = n2
                    ijp = 0
                    js = n2*lda
                    do j = 0,n1 - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,n1
                       do ij = i,i + (n1 + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    ijp = 0
                    jp = 0
                    do j = 0,k - 1
                       do i = j,n - 1
                          ij = 1 + i + jp
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       jp = jp + lda
                    end do
                    do i = 0,k - 1
                       do j = i,k - 1
                          ij = i + j*lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    ijp = 0
                    do j = 0,k - 1
                       ij = k + 1 + j
                       do i = 0,j
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                          ij = ij + lda
                       end do
                    end do
                    js = 0
                    do j = k,n - 1
                       ij = js
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    ijp = 0
                    do i = 0,k - 1
                       do ij = i + (i + 1)*lda, (n + 1)*lda - 1,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                    js = 0
                    do j = 0,k - 1
                       do ij = js,js + k - j - 1
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda + 1
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    ijp = 0
                    js = (k + 1)*lda
                    do j = 0,k - 1
                       do ij = js,js + j
                          arf(ij) = ap(ijp)
                          ijp = ijp + 1
                       end do
                       js = js + lda
                    end do
                    do i = 0,k - 1
                       do ij = i,i + (k + i)*lda,lda
                          arf(ij) = conjg(ap(ijp))
                          ijp = ijp + 1
                       end do
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wtpttf

     !> CTPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_ctpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(sp),intent(out) :: a(lda,*)
           complex(sp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CTPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_ctpttr
     !> ZTPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_ztpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(dp),intent(out) :: a(lda,*)
           complex(dp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZTPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_ztpttr
     !> WTPTTR: copies a triangular matrix A from standard packed format (TP)
     !> to standard full format (TR).

     pure subroutine la_wtpttr(uplo,n,ap,a,lda,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(qp),intent(out) :: a(lda,*)
           complex(qp),intent(in) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WTPTTR',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    a(i,j) = ap(k)
                 end do
              end do
           end if
           return
     end subroutine la_wtpttr

     !> CTRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_ctrttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(sp),intent(in) :: a(0:lda - 1,0:*)
           complex(sp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CTRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    arf(0) = a(0,0)
                 else
                    arf(0) = conjg(a(0,0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = conjg(a(n2 + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = conjg(a(j - n1,l))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda=n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = conjg(a(n2 + j,l))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = conjg(a(k + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = conjg(a(j - k,l))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = conjg(a(k + 1 + j,l))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ctrttf
     !> ZTRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_ztrttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(dp),intent(in) :: a(0:lda - 1,0:*)
           complex(dp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZTRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    arf(0) = a(0,0)
                 else
                    arf(0) = conjg(a(0,0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = conjg(a(n2 + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = conjg(a(j - n1,l))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda=n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = conjg(a(n2 + j,l))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = conjg(a(k + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = conjg(a(j - k,l))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = conjg(a(k + 1 + j,l))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ztrttf
     !> WTRTTF: copies a triangular matrix A from standard full format (TR)
     !> to rectangular full packed format (TF) .

     pure subroutine la_wtrttf(transr,uplo,n,a,lda,arf,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(qp),intent(in) :: a(0:lda - 1,0:*)
           complex(qp),intent(out) :: arf(0:*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: i,ij,j,k,l,n1,n2,nt,nx2,np1x2
           ! Intrinsic Functions
           intrinsic :: conjg,max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WTRTTF',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) then
              if (n == 1) then
                 if (normaltransr) then
                    arf(0) = a(0,0)
                 else
                    arf(0) = conjg(a(0,0))
                 end if
              end if
              return
           end if
           ! size of array arf(1:2,0:nt-1)
           nt = n*(n + 1)/2
           ! set n1 and n2 depending on lower: for n even n1=n2=k
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! if n is odd, set nisodd = .true., lda=n+1 and a is (n+1)--by--k2.
           ! if n is even, set k = n/2 and nisodd = .false., lda=n and a is
           ! n--by--(n+1)/2.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
              if (.not. lower) np1x2 = n + n + 2
           else
              nisodd = .true.
              if (.not. lower) nx2 = n + n
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1); lda=n
                    ij = 0
                    do j = 0,n2
                       do i = n1,n2 + j
                          arf(ij) = conjg(a(n2 + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0); lda=n
                    ij = nt - n
                    do j = n - 1,n1,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - n1,n1 - 1
                          arf(ij) = conjg(a(j - n1,l))
                          ij = ij + 1
                       end do
                       ij = ij - nx2
                    end do
                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0,0) , t2 -> a(1,0) , s -> a(0,n1)
                    ! t1 -> a(0+0) , t2 -> a(1+0) , s -> a(0+n1*n1); lda=n1
                    ij = 0
                    do j = 0,n2 - 1
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = n1 + j,n - 1
                          arf(ij) = a(i,n1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = n2,n - 1
                       do i = 0,n1 - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0,n1+1), t2 -> a(0,n1), s -> a(0,0)
                    ! t1 -> a(n2*n2), t2 -> a(n1*n2), s -> a(0); lda=n2
                    ij = 0
                    do j = 0,n1
                       do i = n1,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,n1 - 1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = n2 + j,n - 1
                          arf(ij) = conjg(a(n2 + j,l))
                          ij = ij + 1
                       end do
                    end do
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1); lda=n+1
                    ij = 0
                    do j = 0,k - 1
                       do i = k,k + j
                          arf(ij) = conjg(a(k + j,i))
                          ij = ij + 1
                       end do
                       do i = j,n - 1
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0); lda=n+1
                    ij = nt - n - 1
                    do j = n - 1,k,-1
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = j - k,k - 1
                          arf(ij) = conjg(a(j - k,l))
                          ij = ij + 1
                       end do
                       ij = ij - np1x2
                    end do
                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,1) , t2 -> a(0,0) , s -> a(0,k+1) :
                    ! t1 -> a(0+k) , t2 -> a(0+0) , s -> a(0+k*(k+1)); lda=k
                    ij = 0
                    j = k
                    do i = k,n - 1
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                       do i = k + 1 + j,n - 1
                          arf(ij) = a(i,k + 1 + j)
                          ij = ij + 1
                       end do
                    end do
                    do j = k - 1,n - 1
                       do i = 0,k - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                 else
                    ! srpa for upper, transpose and n is even (see paper, a=b)
                    ! t1 -> a(0,k+1) , t2 -> a(0,k) , s -> a(0,0)
                    ! t1 -> a(0+k*(k+1)) , t2 -> a(0+k*k) , s -> a(0+0)); lda=k
                    ij = 0
                    do j = 0,k
                       do i = k,n - 1
                          arf(ij) = conjg(a(j,i))
                          ij = ij + 1
                       end do
                    end do
                    do j = 0,k - 2
                       do i = 0,j
                          arf(ij) = a(i,j)
                          ij = ij + 1
                       end do
                       do l = k + 1 + j,n - 1
                          arf(ij) = conjg(a(k + 1 + j,l))
                          ij = ij + 1
                       end do
                    end do
                    ! note that here j = k-1
                    do i = 0,j
                       arf(ij) = a(i,j)
                       ij = ij + 1
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wtrttf

     !> CTRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_ctrttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CTRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_ctrttp
     !> ZTRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_ztrttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZTRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_ztrttp
     !> WTRTTP: copies a triangular matrix A from full format (TR) to standard
     !> packed format (TP).

     pure subroutine la_wtrttp(uplo,n,a,lda,ap,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n,lda
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: ap(*)
        ! =====================================================================
           ! Parameters
           ! Local Scalars
           logical(lk) :: lower
           integer(ilp) :: i,j,k
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WTRTTP',-info)
              return
           end if
           if (lower) then
              k = 0
              do j = 1,n
                 do i = j,n
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           else
              k = 0
              do j = 1,n
                 do i = 1,j
                    k = k + 1
                    ap(k) = a(i,j)
                 end do
              end do
           end if
           return
     end subroutine la_wtrttp

end module la_lapack_blas_like_base
