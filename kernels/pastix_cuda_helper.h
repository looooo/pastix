#ifndef _pastix_cuda_helper_h_
#define _pastix_cuda_helper_h_

#ifdef PASTIX_WITH_CUDA
#  ifdef PRECISION_z
#    define CU_FLOAT           cuDoubleComplex
#    define CU_FLOAT_INIT(r,i) (make_cuDoubleComplex(r,i))
#  endif
#  ifdef PRECISION_c
#    define CU_FLOAT           cuFloatComplex
#    define CU_FLOAT_INIT(r,i) (make_cuFloatComplex(r,i))
#  endif
#  ifdef PRECISION_d
#    define CU_FLOAT double
#    define CU_FLOAT_INIT(r,i) (r)
#  endif
#  ifdef PRECISION_s
#    define CU_FLOAT float
#    define CU_FLOAT_INIT(r,i) (r)
#  endif
#endif /* PASTIX_WITH_CUDA */

#endif /* _pastix_cuda_helper_h_ */
