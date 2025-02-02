!> Matrix and Vector norms
module la_norms
     use la_constants
     use la_blas,only: nrm2
     use la_lapack,only: lange
     use la_state_type
     implicit none(type,external)
     private
     
     public :: norm,get_norm,mnorm

     character(*),parameter :: this = 'norm'
     
     !> List of internal norm flags
     integer(ilp),parameter :: NORM_ONE = 1_ilp
     integer(ilp),parameter :: NORM_TWO = 2_ilp
     integer(ilp),parameter :: NORM_POW_FIRST = 3_ilp
     integer(ilp),parameter :: NORM_INF = +huge(0_ilp) ! infinity norm
     integer(ilp),parameter :: NORM_POW_LAST = NORM_INF - 1_ilp
     integer(ilp),parameter :: NORM_MINUSINF = -huge(0_ilp)
     
     !> List of *LANGE norm flags
     character,parameter :: LANGE_NORM_MAT = 'M' !> maxval(sum(abs(a)))   over whole matrix: unused
     character,parameter :: LANGE_NORM_ONE = '1' !> maxval(sum(abs(a),1)) over columns
     character,parameter :: LANGE_NORM_INF = 'I' !> maxval(sum(abs(a),2)) over rows
     character,parameter :: LANGE_NORM_TWO = 'E' !> "Euclidean" or "Frobenius"
     
     !> @brief Compute the norm of a vector or matrix using LAPACK-based routines.
     !!
     !! Return one of several scalar norm metrics of a real or complex input array A, that can have any rank. 
     !! For generic rank-n arrays, the scalar norm over the whole array is returned by default. 
     !! If \f$ n \geq 2 \f$ and the optional input dimension `dim` is specified, a rank \f$ n-1 \f$ array is returned with dimension `dim` collapsed, 
     !! containing all 1D array norms evaluated along dimension `dim` only.
     !!
     !! Norm type input is mandatory, and it is provided via the `order` argument. This can be provided as either an integer value or a character string. 
     !! Allowed metrics are:
     !! - 1-norm \f$ \sum_i |a_i| \f$: `order = 1` or `order = "1"`
     !! - Euclidean norm \f$ \sqrt{\sum_i a_i^2} \f$: `order = 2` or `order = "2"`
     !! - p-norm \f$ \left( \sum_i |a_i|^p \right)^{1/p} \f$: integer `order >= 3` 
     !! - Infinity norm \f$ \max_i |a_i| \f$: `order = huge(0)` or `"inf"`
     !! - Minus-infinity norm \f$ \min_i |a_i| \f$: `order = -huge(0)` or `"-inf"`
     !!
     !! @param[in] a The input vector or matrix. It may have rank 1 (vector) or higher.
     !! @param[in] order The order of the norm to compute, typically one of the allowed metrics.
     !! @param[in] dim (Optional) The dimension along which to compute the norm, 
     !!               applicable for array norms reducing rank.
     !! @param[out] err (Optional) A state return flag. If an error occurs and `err` is not provided,
     !!                 the function will stop execution.
     !!
     !! @return The computed norm value. If `dim` is specified, the result is a lower-rank array;
     !!         otherwise, it is a scalar.
     !!
     !! @note This interface utilizes LAPACK routines for efficient computation, ensuring numerical stability.
     !! @warning If invalid input values (such as negative norms) are provided, the behavior is undefined.
     !!   
     interface norm
        module procedure la_norm_1D_order_char_s
        module procedure la_norm_1D_order_err_char_s
        module procedure la_norm_2D_order_char_s
        module procedure la_norm_2D_order_err_char_s
        module procedure la_norm_3D_order_char_s
        module procedure la_norm_3D_order_err_char_s
        module procedure la_norm_4D_order_char_s
        module procedure la_norm_4D_order_err_char_s
        module procedure la_norm_5D_order_char_s
        module procedure la_norm_5D_order_err_char_s
        module procedure la_norm_6D_order_char_s
        module procedure la_norm_6D_order_err_char_s
        module procedure la_norm_7D_order_char_s
        module procedure la_norm_7D_order_err_char_s
        module procedure la_norm_2D_to_1D_char_s
        module procedure la_norm_2D_to_1D_err_char_s
        module procedure la_norm_3D_to_2D_char_s
        module procedure la_norm_3D_to_2D_err_char_s
        module procedure la_norm_4D_to_3D_char_s
        module procedure la_norm_4D_to_3D_err_char_s
        module procedure la_norm_5D_to_4D_char_s
        module procedure la_norm_5D_to_4D_err_char_s
        module procedure la_norm_6D_to_5D_char_s
        module procedure la_norm_6D_to_5D_err_char_s
        module procedure la_norm_7D_to_6D_char_s
        module procedure la_norm_7D_to_6D_err_char_s
        module procedure la_norm_1D_order_int_s
        module procedure la_norm_1D_order_err_int_s
        module procedure la_norm_2D_order_int_s
        module procedure la_norm_2D_order_err_int_s
        module procedure la_norm_3D_order_int_s
        module procedure la_norm_3D_order_err_int_s
        module procedure la_norm_4D_order_int_s
        module procedure la_norm_4D_order_err_int_s
        module procedure la_norm_5D_order_int_s
        module procedure la_norm_5D_order_err_int_s
        module procedure la_norm_6D_order_int_s
        module procedure la_norm_6D_order_err_int_s
        module procedure la_norm_7D_order_int_s
        module procedure la_norm_7D_order_err_int_s
        module procedure la_norm_2D_to_1D_int_s
        module procedure la_norm_2D_to_1D_err_int_s
        module procedure la_norm_3D_to_2D_int_s
        module procedure la_norm_3D_to_2D_err_int_s
        module procedure la_norm_4D_to_3D_int_s
        module procedure la_norm_4D_to_3D_err_int_s
        module procedure la_norm_5D_to_4D_int_s
        module procedure la_norm_5D_to_4D_err_int_s
        module procedure la_norm_6D_to_5D_int_s
        module procedure la_norm_6D_to_5D_err_int_s
        module procedure la_norm_7D_to_6D_int_s
        module procedure la_norm_7D_to_6D_err_int_s
        module procedure la_norm_1D_order_char_d
        module procedure la_norm_1D_order_err_char_d
        module procedure la_norm_2D_order_char_d
        module procedure la_norm_2D_order_err_char_d
        module procedure la_norm_3D_order_char_d
        module procedure la_norm_3D_order_err_char_d
        module procedure la_norm_4D_order_char_d
        module procedure la_norm_4D_order_err_char_d
        module procedure la_norm_5D_order_char_d
        module procedure la_norm_5D_order_err_char_d
        module procedure la_norm_6D_order_char_d
        module procedure la_norm_6D_order_err_char_d
        module procedure la_norm_7D_order_char_d
        module procedure la_norm_7D_order_err_char_d
        module procedure la_norm_2D_to_1D_char_d
        module procedure la_norm_2D_to_1D_err_char_d
        module procedure la_norm_3D_to_2D_char_d
        module procedure la_norm_3D_to_2D_err_char_d
        module procedure la_norm_4D_to_3D_char_d
        module procedure la_norm_4D_to_3D_err_char_d
        module procedure la_norm_5D_to_4D_char_d
        module procedure la_norm_5D_to_4D_err_char_d
        module procedure la_norm_6D_to_5D_char_d
        module procedure la_norm_6D_to_5D_err_char_d
        module procedure la_norm_7D_to_6D_char_d
        module procedure la_norm_7D_to_6D_err_char_d
        module procedure la_norm_1D_order_int_d
        module procedure la_norm_1D_order_err_int_d
        module procedure la_norm_2D_order_int_d
        module procedure la_norm_2D_order_err_int_d
        module procedure la_norm_3D_order_int_d
        module procedure la_norm_3D_order_err_int_d
        module procedure la_norm_4D_order_int_d
        module procedure la_norm_4D_order_err_int_d
        module procedure la_norm_5D_order_int_d
        module procedure la_norm_5D_order_err_int_d
        module procedure la_norm_6D_order_int_d
        module procedure la_norm_6D_order_err_int_d
        module procedure la_norm_7D_order_int_d
        module procedure la_norm_7D_order_err_int_d
        module procedure la_norm_2D_to_1D_int_d
        module procedure la_norm_2D_to_1D_err_int_d
        module procedure la_norm_3D_to_2D_int_d
        module procedure la_norm_3D_to_2D_err_int_d
        module procedure la_norm_4D_to_3D_int_d
        module procedure la_norm_4D_to_3D_err_int_d
        module procedure la_norm_5D_to_4D_int_d
        module procedure la_norm_5D_to_4D_err_int_d
        module procedure la_norm_6D_to_5D_int_d
        module procedure la_norm_6D_to_5D_err_int_d
        module procedure la_norm_7D_to_6D_int_d
        module procedure la_norm_7D_to_6D_err_int_d
        module procedure la_norm_1D_order_char_q
        module procedure la_norm_1D_order_err_char_q
        module procedure la_norm_2D_order_char_q
        module procedure la_norm_2D_order_err_char_q
        module procedure la_norm_3D_order_char_q
        module procedure la_norm_3D_order_err_char_q
        module procedure la_norm_4D_order_char_q
        module procedure la_norm_4D_order_err_char_q
        module procedure la_norm_5D_order_char_q
        module procedure la_norm_5D_order_err_char_q
        module procedure la_norm_6D_order_char_q
        module procedure la_norm_6D_order_err_char_q
        module procedure la_norm_7D_order_char_q
        module procedure la_norm_7D_order_err_char_q
        module procedure la_norm_2D_to_1D_char_q
        module procedure la_norm_2D_to_1D_err_char_q
        module procedure la_norm_3D_to_2D_char_q
        module procedure la_norm_3D_to_2D_err_char_q
        module procedure la_norm_4D_to_3D_char_q
        module procedure la_norm_4D_to_3D_err_char_q
        module procedure la_norm_5D_to_4D_char_q
        module procedure la_norm_5D_to_4D_err_char_q
        module procedure la_norm_6D_to_5D_char_q
        module procedure la_norm_6D_to_5D_err_char_q
        module procedure la_norm_7D_to_6D_char_q
        module procedure la_norm_7D_to_6D_err_char_q
        module procedure la_norm_1D_order_int_q
        module procedure la_norm_1D_order_err_int_q
        module procedure la_norm_2D_order_int_q
        module procedure la_norm_2D_order_err_int_q
        module procedure la_norm_3D_order_int_q
        module procedure la_norm_3D_order_err_int_q
        module procedure la_norm_4D_order_int_q
        module procedure la_norm_4D_order_err_int_q
        module procedure la_norm_5D_order_int_q
        module procedure la_norm_5D_order_err_int_q
        module procedure la_norm_6D_order_int_q
        module procedure la_norm_6D_order_err_int_q
        module procedure la_norm_7D_order_int_q
        module procedure la_norm_7D_order_err_int_q
        module procedure la_norm_2D_to_1D_int_q
        module procedure la_norm_2D_to_1D_err_int_q
        module procedure la_norm_3D_to_2D_int_q
        module procedure la_norm_3D_to_2D_err_int_q
        module procedure la_norm_4D_to_3D_int_q
        module procedure la_norm_4D_to_3D_err_int_q
        module procedure la_norm_5D_to_4D_int_q
        module procedure la_norm_5D_to_4D_err_int_q
        module procedure la_norm_6D_to_5D_int_q
        module procedure la_norm_6D_to_5D_err_int_q
        module procedure la_norm_7D_to_6D_int_q
        module procedure la_norm_7D_to_6D_err_int_q
        module procedure la_norm_1D_order_char_c
        module procedure la_norm_1D_order_err_char_c
        module procedure la_norm_2D_order_char_c
        module procedure la_norm_2D_order_err_char_c
        module procedure la_norm_3D_order_char_c
        module procedure la_norm_3D_order_err_char_c
        module procedure la_norm_4D_order_char_c
        module procedure la_norm_4D_order_err_char_c
        module procedure la_norm_5D_order_char_c
        module procedure la_norm_5D_order_err_char_c
        module procedure la_norm_6D_order_char_c
        module procedure la_norm_6D_order_err_char_c
        module procedure la_norm_7D_order_char_c
        module procedure la_norm_7D_order_err_char_c
        !> Array norms: complex(sp)
        module procedure la_norm_2D_to_1D_char_c
        module procedure la_norm_2D_to_1D_err_char_c
        module procedure la_norm_3D_to_2D_char_c
        module procedure la_norm_3D_to_2D_err_char_c
        module procedure la_norm_4D_to_3D_char_c
        module procedure la_norm_4D_to_3D_err_char_c
        module procedure la_norm_5D_to_4D_char_c
        module procedure la_norm_5D_to_4D_err_char_c
        module procedure la_norm_6D_to_5D_char_c
        module procedure la_norm_6D_to_5D_err_char_c
        module procedure la_norm_7D_to_6D_char_c
        module procedure la_norm_7D_to_6D_err_char_c
        !> Scalar norms: complex(sp)
        module procedure la_norm_1D_order_int_c
        module procedure la_norm_1D_order_err_int_c
        module procedure la_norm_2D_order_int_c
        module procedure la_norm_2D_order_err_int_c
        module procedure la_norm_3D_order_int_c
        module procedure la_norm_3D_order_err_int_c
        module procedure la_norm_4D_order_int_c
        module procedure la_norm_4D_order_err_int_c
        module procedure la_norm_5D_order_int_c
        module procedure la_norm_5D_order_err_int_c
        module procedure la_norm_6D_order_int_c
        module procedure la_norm_6D_order_err_int_c
        module procedure la_norm_7D_order_int_c
        module procedure la_norm_7D_order_err_int_c
        !> Array norms: complex(sp)
        module procedure la_norm_2D_to_1D_int_c
        module procedure la_norm_2D_to_1D_err_int_c
        module procedure la_norm_3D_to_2D_int_c
        module procedure la_norm_3D_to_2D_err_int_c
        module procedure la_norm_4D_to_3D_int_c
        module procedure la_norm_4D_to_3D_err_int_c
        module procedure la_norm_5D_to_4D_int_c
        module procedure la_norm_5D_to_4D_err_int_c
        module procedure la_norm_6D_to_5D_int_c
        module procedure la_norm_6D_to_5D_err_int_c
        module procedure la_norm_7D_to_6D_int_c
        module procedure la_norm_7D_to_6D_err_int_c
        !> Scalar norms: complex(dp)
        module procedure la_norm_1D_order_char_z
        module procedure la_norm_1D_order_err_char_z
        module procedure la_norm_2D_order_char_z
        module procedure la_norm_2D_order_err_char_z
        module procedure la_norm_3D_order_char_z
        module procedure la_norm_3D_order_err_char_z
        module procedure la_norm_4D_order_char_z
        module procedure la_norm_4D_order_err_char_z
        module procedure la_norm_5D_order_char_z
        module procedure la_norm_5D_order_err_char_z
        module procedure la_norm_6D_order_char_z
        module procedure la_norm_6D_order_err_char_z
        module procedure la_norm_7D_order_char_z
        module procedure la_norm_7D_order_err_char_z
        !> Array norms: complex(dp)
        module procedure la_norm_2D_to_1D_char_z
        module procedure la_norm_2D_to_1D_err_char_z
        module procedure la_norm_3D_to_2D_char_z
        module procedure la_norm_3D_to_2D_err_char_z
        module procedure la_norm_4D_to_3D_char_z
        module procedure la_norm_4D_to_3D_err_char_z
        module procedure la_norm_5D_to_4D_char_z
        module procedure la_norm_5D_to_4D_err_char_z
        module procedure la_norm_6D_to_5D_char_z
        module procedure la_norm_6D_to_5D_err_char_z
        module procedure la_norm_7D_to_6D_char_z
        module procedure la_norm_7D_to_6D_err_char_z
        !> Scalar norms: complex(dp)
        module procedure la_norm_1D_order_int_z
        module procedure la_norm_1D_order_err_int_z
        module procedure la_norm_2D_order_int_z
        module procedure la_norm_2D_order_err_int_z
        module procedure la_norm_3D_order_int_z
        module procedure la_norm_3D_order_err_int_z
        module procedure la_norm_4D_order_int_z
        module procedure la_norm_4D_order_err_int_z
        module procedure la_norm_5D_order_int_z
        module procedure la_norm_5D_order_err_int_z
        module procedure la_norm_6D_order_int_z
        module procedure la_norm_6D_order_err_int_z
        module procedure la_norm_7D_order_int_z
        module procedure la_norm_7D_order_err_int_z
        !> Array norms: complex(dp)
        module procedure la_norm_2D_to_1D_int_z
        module procedure la_norm_2D_to_1D_err_int_z
        module procedure la_norm_3D_to_2D_int_z
        module procedure la_norm_3D_to_2D_err_int_z
        module procedure la_norm_4D_to_3D_int_z
        module procedure la_norm_4D_to_3D_err_int_z
        module procedure la_norm_5D_to_4D_int_z
        module procedure la_norm_5D_to_4D_err_int_z
        module procedure la_norm_6D_to_5D_int_z
        module procedure la_norm_6D_to_5D_err_int_z
        module procedure la_norm_7D_to_6D_int_z
        module procedure la_norm_7D_to_6D_err_int_z
        !> Scalar norms: complex(qp)
        module procedure la_norm_1D_order_char_w
        module procedure la_norm_1D_order_err_char_w
        module procedure la_norm_2D_order_char_w
        module procedure la_norm_2D_order_err_char_w
        module procedure la_norm_3D_order_char_w
        module procedure la_norm_3D_order_err_char_w
        module procedure la_norm_4D_order_char_w
        module procedure la_norm_4D_order_err_char_w
        module procedure la_norm_5D_order_char_w
        module procedure la_norm_5D_order_err_char_w
        module procedure la_norm_6D_order_char_w
        module procedure la_norm_6D_order_err_char_w
        module procedure la_norm_7D_order_char_w
        module procedure la_norm_7D_order_err_char_w
        !> Array norms: complex(qp)
        module procedure la_norm_2D_to_1D_char_w
        module procedure la_norm_2D_to_1D_err_char_w
        module procedure la_norm_3D_to_2D_char_w
        module procedure la_norm_3D_to_2D_err_char_w
        module procedure la_norm_4D_to_3D_char_w
        module procedure la_norm_4D_to_3D_err_char_w
        module procedure la_norm_5D_to_4D_char_w
        module procedure la_norm_5D_to_4D_err_char_w
        module procedure la_norm_6D_to_5D_char_w
        module procedure la_norm_6D_to_5D_err_char_w
        module procedure la_norm_7D_to_6D_char_w
        module procedure la_norm_7D_to_6D_err_char_w
        !> Scalar norms: complex(qp)
        module procedure la_norm_1D_order_int_w
        module procedure la_norm_1D_order_err_int_w
        module procedure la_norm_2D_order_int_w
        module procedure la_norm_2D_order_err_int_w
        module procedure la_norm_3D_order_int_w
        module procedure la_norm_3D_order_err_int_w
        module procedure la_norm_4D_order_int_w
        module procedure la_norm_4D_order_err_int_w
        module procedure la_norm_5D_order_int_w
        module procedure la_norm_5D_order_err_int_w
        module procedure la_norm_6D_order_int_w
        module procedure la_norm_6D_order_err_int_w
        module procedure la_norm_7D_order_int_w
        module procedure la_norm_7D_order_err_int_w
        module procedure la_norm_2D_to_1D_int_w
        module procedure la_norm_2D_to_1D_err_int_w
        module procedure la_norm_3D_to_2D_int_w
        module procedure la_norm_3D_to_2D_err_int_w
        module procedure la_norm_4D_to_3D_int_w
        module procedure la_norm_4D_to_3D_err_int_w
        module procedure la_norm_5D_to_4D_int_w
        module procedure la_norm_5D_to_4D_err_int_w
        module procedure la_norm_6D_to_5D_int_w
        module procedure la_norm_6D_to_5D_err_int_w
        module procedure la_norm_7D_to_6D_int_w
        module procedure la_norm_7D_to_6D_err_int_w
     end interface norm
     
     !> Vector norm: subroutine interface
     interface get_norm
            !> Scalar norms: real(sp)
            module procedure norm_1D_char_s
            module procedure norm_2D_char_s
            module procedure norm_3D_char_s
            module procedure norm_4D_char_s
            module procedure norm_5D_char_s
            module procedure norm_6D_char_s
            module procedure norm_7D_char_s
            !> Array norms: real(sp)
            module procedure norm_2D_to_1D_char_s
            module procedure norm_3D_to_2D_char_s
            module procedure norm_4D_to_3D_char_s
            module procedure norm_5D_to_4D_char_s
            module procedure norm_6D_to_5D_char_s
            module procedure norm_7D_to_6D_char_s
            !> Scalar norms: real(sp)
            module procedure norm_1D_int_s
            module procedure norm_2D_int_s
            module procedure norm_3D_int_s
            module procedure norm_4D_int_s
            module procedure norm_5D_int_s
            module procedure norm_6D_int_s
            module procedure norm_7D_int_s
            !> Array norms: real(sp)
            module procedure norm_2D_to_1D_int_s
            module procedure norm_3D_to_2D_int_s
            module procedure norm_4D_to_3D_int_s
            module procedure norm_5D_to_4D_int_s
            module procedure norm_6D_to_5D_int_s
            module procedure norm_7D_to_6D_int_s
            !> Scalar norms: real(dp)
            module procedure norm_1D_char_d
            module procedure norm_2D_char_d
            module procedure norm_3D_char_d
            module procedure norm_4D_char_d
            module procedure norm_5D_char_d
            module procedure norm_6D_char_d
            module procedure norm_7D_char_d
            !> Array norms: real(dp)
            module procedure norm_2D_to_1D_char_d
            module procedure norm_3D_to_2D_char_d
            module procedure norm_4D_to_3D_char_d
            module procedure norm_5D_to_4D_char_d
            module procedure norm_6D_to_5D_char_d
            module procedure norm_7D_to_6D_char_d
            !> Scalar norms: real(dp)
            module procedure norm_1D_int_d
            module procedure norm_2D_int_d
            module procedure norm_3D_int_d
            module procedure norm_4D_int_d
            module procedure norm_5D_int_d
            module procedure norm_6D_int_d
            module procedure norm_7D_int_d
            !> Array norms: real(dp)
            module procedure norm_2D_to_1D_int_d
            module procedure norm_3D_to_2D_int_d
            module procedure norm_4D_to_3D_int_d
            module procedure norm_5D_to_4D_int_d
            module procedure norm_6D_to_5D_int_d
            module procedure norm_7D_to_6D_int_d
            !> Scalar norms: real(qp)
            module procedure norm_1D_char_q
            module procedure norm_2D_char_q
            module procedure norm_3D_char_q
            module procedure norm_4D_char_q
            module procedure norm_5D_char_q
            module procedure norm_6D_char_q
            module procedure norm_7D_char_q
            !> Array norms: real(qp)
            module procedure norm_2D_to_1D_char_q
            module procedure norm_3D_to_2D_char_q
            module procedure norm_4D_to_3D_char_q
            module procedure norm_5D_to_4D_char_q
            module procedure norm_6D_to_5D_char_q
            module procedure norm_7D_to_6D_char_q
            !> Scalar norms: real(qp)
            module procedure norm_1D_int_q
            module procedure norm_2D_int_q
            module procedure norm_3D_int_q
            module procedure norm_4D_int_q
            module procedure norm_5D_int_q
            module procedure norm_6D_int_q
            module procedure norm_7D_int_q
            !> Array norms: real(qp)
            module procedure norm_2D_to_1D_int_q
            module procedure norm_3D_to_2D_int_q
            module procedure norm_4D_to_3D_int_q
            module procedure norm_5D_to_4D_int_q
            module procedure norm_6D_to_5D_int_q
            module procedure norm_7D_to_6D_int_q
            !> Scalar norms: complex(sp)
            module procedure norm_1D_char_c
            module procedure norm_2D_char_c
            module procedure norm_3D_char_c
            module procedure norm_4D_char_c
            module procedure norm_5D_char_c
            module procedure norm_6D_char_c
            module procedure norm_7D_char_c
            !> Array norms: complex(sp)
            module procedure norm_2D_to_1D_char_c
            module procedure norm_3D_to_2D_char_c
            module procedure norm_4D_to_3D_char_c
            module procedure norm_5D_to_4D_char_c
            module procedure norm_6D_to_5D_char_c
            module procedure norm_7D_to_6D_char_c
            !> Scalar norms: complex(sp)
            module procedure norm_1D_int_c
            module procedure norm_2D_int_c
            module procedure norm_3D_int_c
            module procedure norm_4D_int_c
            module procedure norm_5D_int_c
            module procedure norm_6D_int_c
            module procedure norm_7D_int_c
            !> Array norms: complex(sp)
            module procedure norm_2D_to_1D_int_c
            module procedure norm_3D_to_2D_int_c
            module procedure norm_4D_to_3D_int_c
            module procedure norm_5D_to_4D_int_c
            module procedure norm_6D_to_5D_int_c
            module procedure norm_7D_to_6D_int_c
            !> Scalar norms: complex(dp)
            module procedure norm_1D_char_z
            module procedure norm_2D_char_z
            module procedure norm_3D_char_z
            module procedure norm_4D_char_z
            module procedure norm_5D_char_z
            module procedure norm_6D_char_z
            module procedure norm_7D_char_z
            !> Array norms: complex(dp)
            module procedure norm_2D_to_1D_char_z
            module procedure norm_3D_to_2D_char_z
            module procedure norm_4D_to_3D_char_z
            module procedure norm_5D_to_4D_char_z
            module procedure norm_6D_to_5D_char_z
            module procedure norm_7D_to_6D_char_z
            !> Scalar norms: complex(dp)
            module procedure norm_1D_int_z
            module procedure norm_2D_int_z
            module procedure norm_3D_int_z
            module procedure norm_4D_int_z
            module procedure norm_5D_int_z
            module procedure norm_6D_int_z
            module procedure norm_7D_int_z
            !> Array norms: complex(dp)
            module procedure norm_2D_to_1D_int_z
            module procedure norm_3D_to_2D_int_z
            module procedure norm_4D_to_3D_int_z
            module procedure norm_5D_to_4D_int_z
            module procedure norm_6D_to_5D_int_z
            module procedure norm_7D_to_6D_int_z
            !> Scalar norms: complex(qp)
            module procedure norm_1D_char_w
            module procedure norm_2D_char_w
            module procedure norm_3D_char_w
            module procedure norm_4D_char_w
            module procedure norm_5D_char_w
            module procedure norm_6D_char_w
            module procedure norm_7D_char_w
            !> Array norms: complex(qp)
            module procedure norm_2D_to_1D_char_w
            module procedure norm_3D_to_2D_char_w
            module procedure norm_4D_to_3D_char_w
            module procedure norm_5D_to_4D_char_w
            module procedure norm_6D_to_5D_char_w
            module procedure norm_7D_to_6D_char_w
            !> Scalar norms: complex(qp)
            module procedure norm_1D_int_w
            module procedure norm_2D_int_w
            module procedure norm_3D_int_w
            module procedure norm_4D_int_w
            module procedure norm_5D_int_w
            module procedure norm_6D_int_w
            module procedure norm_7D_int_w
            !> Array norms: complex(qp)
            module procedure norm_2D_to_1D_int_w
            module procedure norm_3D_to_2D_int_w
            module procedure norm_4D_to_3D_int_w
            module procedure norm_5D_to_4D_int_w
            module procedure norm_6D_to_5D_int_w
            module procedure norm_7D_to_6D_int_w
     end interface get_norm
     
     !> Matrix norm: function interface
     interface mnorm
         module procedure matrix_norm_char_s
         module procedure matrix_norm_3D_to_1D_char_s
         module procedure matrix_norm_4D_to_2D_char_s
         module procedure matrix_norm_5D_to_3D_char_s
         module procedure matrix_norm_6D_to_4D_char_s
         module procedure matrix_norm_7D_to_5D_char_s
         module procedure matrix_norm_int_s
         module procedure matrix_norm_3D_to_1D_int_s
         module procedure matrix_norm_4D_to_2D_int_s
         module procedure matrix_norm_5D_to_3D_int_s
         module procedure matrix_norm_6D_to_4D_int_s
         module procedure matrix_norm_7D_to_5D_int_s
         module procedure matrix_norm_char_d
         module procedure matrix_norm_3D_to_1D_char_d
         module procedure matrix_norm_4D_to_2D_char_d
         module procedure matrix_norm_5D_to_3D_char_d
         module procedure matrix_norm_6D_to_4D_char_d
         module procedure matrix_norm_7D_to_5D_char_d
         module procedure matrix_norm_int_d
         module procedure matrix_norm_3D_to_1D_int_d
         module procedure matrix_norm_4D_to_2D_int_d
         module procedure matrix_norm_5D_to_3D_int_d
         module procedure matrix_norm_6D_to_4D_int_d
         module procedure matrix_norm_7D_to_5D_int_d
         module procedure matrix_norm_char_q
         module procedure matrix_norm_3D_to_1D_char_q
         module procedure matrix_norm_4D_to_2D_char_q
         module procedure matrix_norm_5D_to_3D_char_q
         module procedure matrix_norm_6D_to_4D_char_q
         module procedure matrix_norm_7D_to_5D_char_q
         module procedure matrix_norm_int_q
         module procedure matrix_norm_3D_to_1D_int_q
         module procedure matrix_norm_4D_to_2D_int_q
         module procedure matrix_norm_5D_to_3D_int_q
         module procedure matrix_norm_6D_to_4D_int_q
         module procedure matrix_norm_7D_to_5D_int_q
         module procedure matrix_norm_char_c
         module procedure matrix_norm_3D_to_1D_char_c
         module procedure matrix_norm_4D_to_2D_char_c
         module procedure matrix_norm_5D_to_3D_char_c
         module procedure matrix_norm_6D_to_4D_char_c
         module procedure matrix_norm_7D_to_5D_char_c
         module procedure matrix_norm_int_c
         module procedure matrix_norm_3D_to_1D_int_c
         module procedure matrix_norm_4D_to_2D_int_c
         module procedure matrix_norm_5D_to_3D_int_c
         module procedure matrix_norm_6D_to_4D_int_c
         module procedure matrix_norm_7D_to_5D_int_c
         module procedure matrix_norm_char_z
         module procedure matrix_norm_3D_to_1D_char_z
         module procedure matrix_norm_4D_to_2D_char_z
         module procedure matrix_norm_5D_to_3D_char_z
         module procedure matrix_norm_6D_to_4D_char_z
         module procedure matrix_norm_7D_to_5D_char_z
         module procedure matrix_norm_int_z
         module procedure matrix_norm_3D_to_1D_int_z
         module procedure matrix_norm_4D_to_2D_int_z
         module procedure matrix_norm_5D_to_3D_int_z
         module procedure matrix_norm_6D_to_4D_int_z
         module procedure matrix_norm_7D_to_5D_int_z
         module procedure matrix_norm_char_w
         module procedure matrix_norm_3D_to_1D_char_w
         module procedure matrix_norm_4D_to_2D_char_w
         module procedure matrix_norm_5D_to_3D_char_w
         module procedure matrix_norm_6D_to_4D_char_w
         module procedure matrix_norm_7D_to_5D_char_w
         module procedure matrix_norm_int_w
         module procedure matrix_norm_3D_to_1D_int_w
         module procedure matrix_norm_4D_to_2D_int_w
         module procedure matrix_norm_5D_to_3D_int_w
         module procedure matrix_norm_6D_to_4D_int_w
         module procedure matrix_norm_7D_to_5D_int_w
     end interface mnorm
     
     interface parse_norm_type
        module procedure parse_norm_type_integer
        module procedure parse_norm_type_character
     end interface parse_norm_type
     
     contains
     
     !> Parse norm type from an integer user input
     pure subroutine parse_norm_type_integer(order,norm_type,err)
        !> User input value
        integer(ilp),intent(in) :: order
        !> Return value: norm type
        integer(ilp),intent(out) :: norm_type
        !> State return flag
        type(la_state),intent(out) :: err
        
        select case (order)
           case (1_ilp)
               norm_type = NORM_ONE
           case (2_ilp)
               norm_type = NORM_TWO
           case (3_ilp:huge(0_ilp) - 1_ilp)
               norm_type = order
           case (huge(0_ilp):)
               norm_type = NORM_INF
           case (:-huge(0_ilp))
               norm_type = NORM_MINUSINF
           
           case default
               norm_type = NORM_ONE
               err = la_state(this,LINALG_ERROR,'Input norm type ',order,' is not recognized.')
        end select
        
     end subroutine parse_norm_type_integer

     pure subroutine parse_norm_type_character(order,norm_type,err)
        !> User input value
        character(len=*),intent(in) :: order
        !> Return value: norm type
        integer(ilp),intent(out) :: norm_type
        !> State return flag
        type(la_state),intent(out) :: err
        
        integer(ilp) :: int_order,read_err
        
        select case (order)
           case ('inf','Inf','INF')
              norm_type = NORM_INF
           case ('-inf','-Inf','-INF')
              norm_type = NORM_MINUSINF
           case ('Euclidean','euclidean','EUCLIDEAN')
              norm_type = NORM_TWO
           case default
            
              ! Check if this input can be read as an integer
              read (order,*,iostat=read_err) int_order
              if (read_err /= 0) then
                 ! Cannot read as an integer
                 norm_type = NORM_ONE
                 err = la_state(this,LINALG_ERROR,'Input norm type ',order,' is not recognized.')
              else
                 call parse_norm_type_integer(int_order,norm_type,err)
              end if

        end select
        
     end subroutine parse_norm_type_character

     !> From a user norm request, generate a *LANGE task command
     pure subroutine lange_task_request(norm_type,lange_task,err)
        !> Parsed matrix norm type
        integer(ilp),intent(in) :: norm_type
        !> LANGE task
        character,intent(out) :: lange_task
        !> Error flag
        type(la_state),intent(inout) :: err
        
        select case (norm_type)
           case (NORM_INF)
              lange_task = LANGE_NORM_INF
           case (NORM_ONE)
              lange_task = LANGE_NORM_ONE
           case (NORM_TWO)
              lange_task = LANGE_NORM_TWO
           case default
              err = la_state(this,LINALG_VALUE_ERROR,'Order ',norm_type,' is not a valid matrix norm input.')
        end select
     end subroutine lange_task_request
               
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_s(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_1D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_s
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_s(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_1D_char_s(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_s(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_2D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_s
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_s(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_2D_char_s(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_s(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_3D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_s
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_s(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_3D_char_s(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_s(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_4D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_s
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_s(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_4D_char_s(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_s(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_5D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_s
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_s(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_5D_char_s(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_s(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_6D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_s
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_s(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_6D_char_s(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_s(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_7D_char_s(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_s
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_s(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_char_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_s
    
    ! Internal implementation
    pure subroutine norm_7D_char_s(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_s

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_s(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_s(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_s(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_s(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_s(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_s(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_s(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_s
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_s

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_s(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(sp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_s
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_s(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_1D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_s
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_s(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_1D_int_s(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(sp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_s(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_2D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_s
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_s(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_2D_int_s(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(sp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_s(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_3D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_s
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_s(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_3D_int_s(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(sp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_s(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_4D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_s
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_s(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_4D_int_s(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_s(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_5D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_s
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_s(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_5D_int_s(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_s(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_6D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_s
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_s(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_6D_int_s(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_s

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_s(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_7D_int_s(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_s
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_s(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_int_s(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_s
    
    ! Internal implementation
    pure subroutine norm_7D_int_s(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_s

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_s(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_s(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_s(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_s(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_s(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_s

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_s(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_s(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_s

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_s(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_s
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_s(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_s

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_s(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(sp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_s
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_s(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        real(sp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_s
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_d(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_1D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_d
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_d(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_1D_char_d(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_d(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_2D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_d
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_d(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_2D_char_d(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_d(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_3D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_d
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_d(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_3D_char_d(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_d(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_4D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_d
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_d(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_4D_char_d(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_d(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_5D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_d
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_d(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_5D_char_d(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_d(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_6D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_d
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_d(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_6D_char_d(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_d(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_7D_char_d(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_d
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_d(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_char_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_d
    
    ! Internal implementation
    pure subroutine norm_7D_char_d(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_d

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_d(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_d(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_d(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_d(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_d(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_d(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_d(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_d
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_d

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_d(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(dp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_d
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_d(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_1D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_d
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_d(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_1D_int_d(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(dp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_d(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_2D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_d
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_d(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_2D_int_d(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(dp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_d(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_3D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_d
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_d(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_3D_int_d(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(dp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_d(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_4D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_d
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_d(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_4D_int_d(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_d(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_5D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_d
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_d(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_5D_int_d(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_d(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_6D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_d
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_d(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_6D_int_d(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_d

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_d(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_7D_int_d(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_d
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_d(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_int_d(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_d
    
    ! Internal implementation
    pure subroutine norm_7D_int_d(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_d

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_d(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_d(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_d(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_d(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_d(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_d

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_d(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_d(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_d

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_d(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_d
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_d(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_d

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_d(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(dp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_d
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_d(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        real(dp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_d
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_q(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_1D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_q
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_q(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_1D_char_q(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_q(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_2D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_q
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_q(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_2D_char_q(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_q(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_3D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_q
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_q(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_3D_char_q(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_q(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_4D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_q
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_q(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_4D_char_q(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_q(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_5D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_q
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_q(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_5D_char_q(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_q(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_6D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_q
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_q(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_6D_char_q(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_q(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_7D_char_q(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_q
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_q(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_char_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_q
    
    ! Internal implementation
    pure subroutine norm_7D_char_q(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_q

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_q(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_q(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_q(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_q(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_q(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_q(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_q(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_q
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_q

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_q(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(qp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_q
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_q(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_1D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_q
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_q(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_1D_int_q(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        real(qp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_q(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_2D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_q
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_q(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_2D_int_q(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        real(qp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_q(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_3D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_q
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_q(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_3D_int_q(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        real(qp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_q(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_4D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_q
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_q(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_4D_int_q(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_q(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_5D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_q
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_q(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_5D_int_q(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_q(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_6D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_q
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_q(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_6D_int_q(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_q

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_q(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_7D_int_q(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_q
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_q(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_int_q(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_q
    
    ! Internal implementation
    pure subroutine norm_7D_int_q(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        real(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = norm2(a)
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_q

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_q(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_q(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_q(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_q(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_q(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_q

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_q(a,order,dim) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_q(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_q

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_q(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_q
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_q(a,nrm,order,dim,err)
        !> Input matrix a[..]
        real(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = norm2(a,dim=dim)
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_q

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_q(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(qp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_q
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_q(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        real(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        real(qp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_q
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_c(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_1D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_c
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_c(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_1D_char_c(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_c(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_2D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_c
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_c(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_2D_char_c(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_c(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_3D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_c
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_c(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_3D_char_c(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_c(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_4D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_c
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_c(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_4D_char_c(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_c(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_5D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_c
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_c(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_5D_char_c(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_c(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_6D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_c
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_c(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_6D_char_c(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_c(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_7D_char_c(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_c
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_c(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_char_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_c
    
    ! Internal implementation
    pure subroutine norm_7D_char_c(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_c

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_c(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_c(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_c(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_c(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_c(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_c(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_c(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_c
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_c

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_c(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(sp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_c
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_c(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_1D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_c
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_c(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_1D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_1D_int_c(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(sp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_c(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_2D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_c
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_c(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_2D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_2D_int_c(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(sp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_c(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_3D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_c
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_c(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_3D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_3D_int_c(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(sp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_c(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_4D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_c
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_c(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_4D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_4D_int_c(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_c(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_5D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_c
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_c(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_5D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_5D_int_c(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_c(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_6D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_c
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_c(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_6D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_6D_int_c(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_c

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_c(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(sp) :: nrm
                                    
        call norm_7D_int_c(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_c
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_c(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm
                
        call norm_7D_int_c(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_c
    
    ! Internal implementation
    pure subroutine norm_7D_int_c(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(sp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),sp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_c

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_c(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_c(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_c(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_c(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_c(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_c

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_c(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_c(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_c

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(sp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_c(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_c
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_c(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(sp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(sp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(sp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),sp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_sp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_c

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_c(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(sp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_sp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_c
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_c(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(sp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(sp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(sp),target :: work1(1)
        real(sp),pointer :: work(:)
        complex(sp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_c
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_z(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_1D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_z
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_z(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_1D_char_z(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_z(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_2D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_z
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_z(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_2D_char_z(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_z(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_3D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_z
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_z(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_3D_char_z(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_z(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_4D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_z
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_z(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_4D_char_z(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_z(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_5D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_z
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_z(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_5D_char_z(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_z(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_6D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_z
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_z(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_6D_char_z(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_z(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_7D_char_z(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_z
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_z(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_char_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_z
    
    ! Internal implementation
    pure subroutine norm_7D_char_z(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_z

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_z(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_z(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_z(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_z(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_z(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_z(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_z(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_z
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_z

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_z(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(dp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_z
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_z(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_1D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_z
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_z(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_1D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_1D_int_z(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(dp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_z(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_2D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_z
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_z(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_2D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_2D_int_z(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(dp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_z(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_3D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_z
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_z(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_3D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_3D_int_z(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(dp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_z(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_4D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_z
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_z(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_4D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_4D_int_z(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_z(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_5D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_z
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_z(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_5D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_5D_int_z(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_z(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_6D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_z
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_z(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_6D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_6D_int_z(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_z

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_z(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(dp) :: nrm
                                    
        call norm_7D_int_z(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_z
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_z(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm
                
        call norm_7D_int_z(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_z
    
    ! Internal implementation
    pure subroutine norm_7D_int_z(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(dp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),dp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_z

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_z(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_z(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_z(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_z(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_z(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_z

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_z(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_z(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_z

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(dp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_z(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_z
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_z(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(dp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(dp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(dp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),dp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_dp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_z

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_z(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(dp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_dp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_z
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_z(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(dp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(dp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(dp),target :: work1(1)
        real(dp),pointer :: work(:)
        complex(dp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_z
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_char_w(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_1D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_char_w
    
    ! Function interface with output error
    function la_norm_1D_order_err_char_w(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_1D_char_w(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_char_w(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_2D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_char_w
    
    ! Function interface with output error
    function la_norm_2D_order_err_char_w(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_2D_char_w(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_char_w(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_3D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_char_w
    
    ! Function interface with output error
    function la_norm_3D_order_err_char_w(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_3D_char_w(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_char_w(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_4D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_char_w
    
    ! Function interface with output error
    function la_norm_4D_order_err_char_w(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_4D_char_w(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_char_w(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_5D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_char_w
    
    ! Function interface with output error
    function la_norm_5D_order_err_char_w(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_5D_char_w(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_char_w(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_6D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_char_w
    
    ! Function interface with output error
    function la_norm_6D_order_err_char_w(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_6D_char_w(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_char_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_char_w(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_7D_char_w(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_char_w
    
    ! Function interface with output error
    function la_norm_7D_order_err_char_w(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_char_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_char_w
    
    ! Internal implementation
    pure subroutine norm_7D_char_w(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_char_w

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and char input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_w(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_char_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_w(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_char_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_w(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_char_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_w(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_char_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_w(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_char_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_char_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_w(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_char_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_char_w(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_char_w
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_char_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_char_w

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_char_w(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(qp) :: nrm
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_char_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_char_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_char_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_char_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_char_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_char_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        character(len=*),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_char_w
    
    !==============================================
    ! Norms : any rank to scalar
    !==============================================

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_1D_order_int_w(a,order) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_1D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_1D_order_int_w
    
    ! Function interface with output error
    function la_norm_1D_order_err_int_w(a,order,err) result(nrm)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_1D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_1D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_1D_int_w(a,nrm,order,err)
        !> Input 1-d matrix a(:)
        complex(qp),intent(in) :: a(:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_1D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_2D_order_int_w(a,order) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_2D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_2D_order_int_w
    
    ! Function interface with output error
    function la_norm_2D_order_err_int_w(a,order,err) result(nrm)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_2D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_2D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_2D_int_w(a,nrm,order,err)
        !> Input 2-d matrix a(:,:)
        complex(qp),intent(in) :: a(:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_3D_order_int_w(a,order) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_3D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_3D_order_int_w
    
    ! Function interface with output error
    function la_norm_3D_order_err_int_w(a,order,err) result(nrm)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_3D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_3D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_3D_int_w(a,nrm,order,err)
        !> Input 3-d matrix a(:,:,:)
        complex(qp),intent(in) :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_4D_order_int_w(a,order) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_4D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_4D_order_int_w
    
    ! Function interface with output error
    function la_norm_4D_order_err_int_w(a,order,err) result(nrm)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_4D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_4D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_4D_int_w(a,nrm,order,err)
        !> Input 4-d matrix a(:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_5D_order_int_w(a,order) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_5D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_5D_order_int_w
    
    ! Function interface with output error
    function la_norm_5D_order_err_int_w(a,order,err) result(nrm)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_5D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_5D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_5D_int_w(a,nrm,order,err)
        !> Input 5-d matrix a(:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_6D_order_int_w(a,order) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_6D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_6D_order_int_w
    
    ! Function interface with output error
    function la_norm_6D_order_err_int_w(a,order,err) result(nrm)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_6D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_6D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_6D_int_w(a,nrm,order,err)
        !> Input 6-d matrix a(:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_int_w

    ! Pure function interface, with order specification. On error, the code will stop
    pure function la_norm_7D_order_int_w(a,order) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Norm of the matrix.
        real(qp) :: nrm
                                    
        call norm_7D_int_w(a,nrm=nrm,order=order)
        
    end function la_norm_7D_order_int_w
    
    ! Function interface with output error
    function la_norm_7D_order_err_int_w(a,order,err) result(nrm)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm
                
        call norm_7D_int_w(a,nrm=nrm,order=order,err=err)
        
    end function la_norm_7D_order_err_int_w
    
    ! Internal implementation
    pure subroutine norm_7D_int_w(a,nrm,order,err)
        !> Input 7-d matrix a(:,:,:,:,:,:,:)
        complex(qp),intent(in) :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        
        intrinsic :: abs,sum,sqrt,norm2,maxval,minval,conjg
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        ! Check matrix size
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a))
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a)),qp))
            case (NORM_INF)
                nrm = maxval(abs(a))
            case (-NORM_INF)
                nrm = minval(abs(a))
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_int_w

    !====================================================================
    ! Norms : any rank to rank-1, with DIM specifier and int input
    !====================================================================

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_2D_to_1D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_w(a,nrm,order,dim)
            
    end function la_norm_2D_to_1D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_2D_to_1D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        
        call norm_2D_to_1D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_2D_to_1D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_2D_to_1D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 2) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_2D_to_1D_int_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_3D_to_2D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_w(a,nrm,order,dim)
            
    end function la_norm_3D_to_2D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_3D_to_2D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        
        call norm_3D_to_2D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_3D_to_2D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_3D_to_2D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 3) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_3D_to_2D_int_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_4D_to_3D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_w(a,nrm,order,dim)
            
    end function la_norm_4D_to_3D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_4D_to_3D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim))
        
        call norm_4D_to_3D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_4D_to_3D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_4D_to_3D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 4) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_4D_to_3D_int_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_5D_to_4D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_w(a,nrm,order,dim)
            
    end function la_norm_5D_to_4D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_5D_to_4D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        
        call norm_5D_to_4D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_5D_to_4D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_5D_to_4D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 5) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_5D_to_4D_int_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_6D_to_5D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_w(a,nrm,order,dim)
            
    end function la_norm_6D_to_5D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_6D_to_5D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim))
        
        call norm_6D_to_5D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_6D_to_5D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_6D_to_5D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 6) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_6D_to_5D_int_w

    ! Pure function interface with DIM specifier. On error, the code will stop
    pure function la_norm_7D_to_6D_int_w(a,order,dim) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_w(a,nrm,order,dim)
            
    end function la_norm_7D_to_6D_int_w

    ! Function interface with DIM specifier and output error state.
    function la_norm_7D_to_6D_err_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        integer(ilp),intent(in) :: dim
        !> Output state return flag.
        type(la_state),intent(out) :: err
        !> Norm of the matrix.
        real(qp) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),merge(size(a,3),&
            & size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),mask=5 < dim),&
            & merge(size(a,6),size(a,7),mask=6 < dim))
        
        call norm_7D_to_6D_int_w(a,nrm,order,dim,err)
            
    end function la_norm_7D_to_6D_err_int_w
    
    ! Internal implementation
    pure subroutine norm_7D_to_6D_int_w(a,nrm,order,dim,err)
        !> Input matrix a[..]
        complex(qp),intent(in),target :: a(:,:,:,:,:,:,:)
        !> Dimension to collapse by computing the norm w.r.t other dimensions
        !  (dim must be defined before it is used for `nrm`)
        integer(ilp),intent(in) :: dim
        !> Norm of the matrix.
        real(qp),intent(out) :: nrm(merge(size(a,1),size(a,2),mask=1 < dim),merge(size(a,2),size(a,3),mask=2 < dim),&
            & merge(size(a,3),size(a,4),mask=3 < dim),merge(size(a,4),size(a,5),mask=4 < dim),merge(size(a,5),size(a,6),&
            & mask=5 < dim),merge(size(a,6),size(a,7),mask=6 < dim))
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: sze,norm_request
        real(qp) :: rorder
        intrinsic :: sum,abs,sqrt,conjg,norm2,maxval,minval
        
        sze = size(a,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (sze <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if

        ! Check dimension choice
        if (dim < 1 .or. dim > 7) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'dimension ',dim, &
                                'is out of rank for shape(a)=',shape(a,kind=ilp))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request
        call parse_norm_type(order,norm_request,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        select case (norm_request)
            case (NORM_ONE)
                nrm = sum(abs(a),dim=dim)
            case (NORM_TWO)
                nrm = sqrt(real(sum(a*conjg(a),dim=dim),qp))
            case (NORM_INF)
                nrm = maxval(abs(a),dim=dim)
            case (-NORM_INF)
                nrm = minval(abs(a),dim=dim)
            case (NORM_POW_FIRST:NORM_POW_LAST)
                rorder = 1.0_qp/norm_request
                nrm = sum(abs(a)**norm_request,dim=dim)**rorder
            case default
                err0 = la_state(this,LINALG_INTERNAL_ERROR,'invalid norm type after checking')
                call err0%handle(err)
        end select
        
    end subroutine norm_7D_to_6D_int_w

    !====================================================================
    ! Matrix norms
    !====================================================================
    
    ! Internal implementation
    function matrix_norm_int_w(a,order,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),target :: a(:,:)
        !> Norm of the matrix.
        real(qp) :: nrm
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: m,n,norm_request
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        
        m = size(a,dim=1,kind=ilp)
        n = size(a,dim=2,kind=ilp)
        
        ! Initialize norm to zero
        nrm = 0.0_qp
        
        if (m <= 0 .or. n <= 0) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'invalid matrix shape: a=', [m,n])
            call err0%handle(err)
            return
        end if

        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! LAPACK interface
        nrm = lange(lange_task,m,n,a,m,work)
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        
    end function matrix_norm_int_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_3D_to_1D_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(3) :: s,spack,perm
        integer(ilp),dimension(3),parameter :: dim_range = [(m,m=1_ilp,3_ilp)]
        integer(ilp) :: j3
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 3)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',3,' matrix norm has invalid dim=',dims)
            allocate (nrm(0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3) = &
        lange(lange_task,m,n,apack(:,:,j3),lda,work)
        end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_3D_to_1D_int_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_4D_to_2D_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(4) :: s,spack,perm
        integer(ilp),dimension(4),parameter :: dim_range = [(m,m=1_ilp,4_ilp)]
        integer(ilp) :: j3,j4
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 4)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',4,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4) = &
        lange(lange_task,m,n,apack(:,:,j3,j4),lda,work)
        end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_4D_to_2D_int_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_5D_to_3D_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(5) :: s,spack,perm
        integer(ilp),dimension(5),parameter :: dim_range = [(m,m=1_ilp,5_ilp)]
        integer(ilp) :: j3,j4,j5
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 5)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',5,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5),lda,work)
        end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_5D_to_3D_int_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_6D_to_4D_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(6) :: s,spack,perm
        integer(ilp),dimension(6),parameter :: dim_range = [(m,m=1_ilp,6_ilp)]
        integer(ilp) :: j3,j4,j5,j6
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 6)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',6,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6),lda,work)
        end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_6D_to_4D_int_w
    
    ! Pure function interface with DIM specifier. On error, the code will stop
    function matrix_norm_7D_to_5D_int_w(a,order,dim,err) result(nrm)
        !> Input matrix a(m,n)
        complex(qp),intent(in),contiguous,target :: a(:,:,:,:,:,:,:)
        !> Norm of the matrix.
        real(qp),allocatable :: nrm(:,:,:,:,:)
        !> Order of the matrix norm being computed.
        integer(ilp),intent(in) :: order
        !> [optional] dimensions of the sub-matrices the norms should be evaluated at (default = [1,2])
        integer(ilp),optional,intent(in) :: dim(2)
        !> [optional] state return flag. On error if not requested, the code will stop
        type(la_state),intent(out),optional :: err
        
        type(la_state) :: err0
        integer(ilp) :: j,m,n,lda,dims(2),norm_request
        integer(ilp),dimension(7) :: s,spack,perm
        integer(ilp),dimension(7),parameter :: dim_range = [(m,m=1_ilp,7_ilp)]
        integer(ilp) :: j3,j4,j5,j6,j7
        logical :: contiguous_data
        character :: lange_task
        real(qp),target :: work1(1)
        real(qp),pointer :: work(:)
        complex(qp),pointer :: apack(:,:,:,:,:,:,:)
        
        ! Get dimensions
        if (present(dim)) then
           dims = dim
        else
           dims = [1,2]
        end if
        
        nullify (apack)

        if (dims(1) == dims(2) .or. .not. all(dims > 0 .and. dims <= 7)) then
            err0 = la_state(this,LINALG_VALUE_ERROR,'Rank-',7,' matrix norm has invalid dim=',dims)
            allocate (nrm(0,0,0,0,0))
            call err0%handle(err)
            return
        end if
        
        ! Check norm request: user + *LANGE support
        call parse_norm_type(order,norm_request,err0)
        call lange_task_request(norm_request,lange_task,err0)
        if (err0%error()) then
            call err0%handle(err)
            return
        end if
        
        ! Input matrix properties
        s = shape(a,kind=ilp)
        
        ! Check if input column data is contiguous
        contiguous_data = dims(1) == 1
        
        ! Matrix norm size
        m = s(dims(1))
        n = s(dims(2))

        ! Get packed data with norm dimensions as 1:2
        if (contiguous_data) then
            
            ! Collapse everything before the 1st dimension as apack's dim #1
            ! Set size==1 for all unused trailing dimensions
            spack = [product(s(:dims(2) - 1)),s(dims(2):), (1_ilp,j=1,dims(2) - 2)]
            
            ! Reshape without moving data
            apack(1:spack(1),1:spack(2),1:spack(3),1:spack(4),1:spack(5),1:spack(6),1:spack(7)) => a
            
        else
            
            ! Dimension permutations to map dims(1),dims(2) => 1:2
            perm = [dims,pack(dim_range,dim_range /= dims(1) .and. dim_range /= dims(2))]
            spack = s(perm)
            apack = reshape(a,shape=spack,order=perm)
            
        end if
            
        if (lange_task == LANGE_NORM_INF) then
            allocate (work(m))
        else
            work => work1
        end if
        
        ! Allocate norm
        allocate (nrm(size(apack,3),size(apack,4),size(apack,5),size(apack,6),size(apack,7)))
        
        lda = size(apack,dim=1,kind=ilp)
        
        ! LAPACK interface
        do j7 = lbound(apack,7),ubound(apack,7)
        do j6 = lbound(apack,6),ubound(apack,6)
        do j5 = lbound(apack,5),ubound(apack,5)
        do j4 = lbound(apack,4),ubound(apack,4)
        do j3 = lbound(apack,3),ubound(apack,3)
        nrm(j3,j4,j5,j6,j7) = &
        lange(lange_task,m,n,apack(:,:,j3,j4,j5,j6,j7),lda,work)
        end do; end do; end do; end do; end do
        
        if (lange_task == LANGE_NORM_INF) deallocate (work)
        if (.not. contiguous_data) deallocate (apack)
        
    end function matrix_norm_7D_to_5D_int_w
    
end module la_norms
