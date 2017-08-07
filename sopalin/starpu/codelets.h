/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelets.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/

#ifndef _CODELETS_H_
#define _CODELETS_H_

#include "common.h"

#ifdef STARPU_CUDA_ASYNC
#define CODELET_CUDA_FLAGS(flags) .cuda_flags = {(flags)},
#else
#define CODELET_CUDA_FLAGS(flags)
#endif

#define CODELETS_ALL( _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ ) \
    struct starpu_codelet cl_##_name_ = {                               \
        .where     = (_original_location_),                             \
        .cpu_funcs[0] = (_cpu_func_name_),                              \
        CODELET_CUDA_FLAGS(_cuda_flags_)                                \
        .cuda_funcs[0] = (_cuda_func_name_),                            \
        .nbuffers  = (_nbuffers_),                                      \
        .model     = &starpu_##_name_##_model,                          \
        .name      = #_name_                                            \
    };

#if defined(PASTIX_STARPU_SIMULATION)

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( _name_, _nbuffers_, (starpu_cpu_func_t) 1, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( _name_, _nbuffers_, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#else

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( _name_, _nbuffers_, cl_##_name_##_cpu, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( _name_, _nbuffers_, cl_##_name_##_cpu, cl_##_name_##_gpu, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#endif

#if !defined(PASTIX_WITH_CUDA)
#undef CODELETS_GPU
#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)       \
    CODELETS_CPU( _name_, _nbuffers_ )
#endif


#endif /* _CODELETS_H_ */
