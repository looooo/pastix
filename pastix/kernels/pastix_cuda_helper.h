#ifndef PASTIX_CUDA_HELPER_H
#define PASTIX_CUDA_HELPER_H

#ifdef PREC_DOUBLE
#  ifdef TYPE_COMPLEX
#    define PRECISION_z
#  else
#    define PRECISION_d
#  endif
#else
#  ifdef TYPE_COMPLEX
#    define PRECISION_c
#  else
#    define PRECISION_s
#  endif
#endif

#ifdef PASTIX_WITH_STARPU
#  ifdef PASTIX_WITH_CUDA
#    ifdef PRECISION_z
#      define CU_FLOAT           cuDoubleComplex
#      define CU_FLOAT_INIT(r,i) (make_cuDoubleComplex(r,i))
#    endif
#    ifdef PRECISION_c
#      define CU_FLOAT           cuFloatComplex
#      define CU_FLOAT_INIT(r,i) (make_cuFloatComplex(r,i))
#    endif
#    ifdef PRECISION_d
#      define CU_FLOAT double
#      define CU_FLOAT_INIT(r,i) (r)
#    endif
#    ifdef PRECISION_s
#      define CU_FLOAT float
#      define CU_FLOAT_INIT(r,i) (r)
#    endif
#  endif /* not FORCE_NO_CUDA */
#endif /* WITH_STARPU */

#endif /* PASTIX_CUDA_HELPER_H */
