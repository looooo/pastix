/**
 *
 * @file codelets.h
 *
 * @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2017-06-01
 *
 **/

#ifndef _codelets_h_
#define _codelets_h_

#include "common.h"

#ifdef STARPU_CUDA_ASYNC
#define CODELET_CUDA_FLAGS(flags) .cuda_flags = {(flags)},
#else
#define CODELET_CUDA_FLAGS(flags)
#endif

#define CODELETS_ALL( _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_) \
    struct starpu_codelet cl_##_name_ = {                               \
        .where     = (_original_location_),                             \
        .cpu_funcs[0] = (_cpu_func_name_),                              \
        CODELET_CUDA_FLAGS(_cuda_flags_)                                \
        .cuda_funcs[0] = (_cuda_func_name_),                            \
        .nbuffers  = (_nbuffers_),                                      \
        .name      = #_name_                                            \
    };

#define CODELETS_ALL_MODEL( _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ ,_perfmodel_) \
    struct starpu_codelet cl_##_name_ = {                               \
        .where     = (_original_location_),                             \
        .cpu_funcs[0] = (_cpu_func_name_),                              \
        CODELET_CUDA_FLAGS(_cuda_flags_)                                \
        .cuda_funcs[0] = (_cuda_func_name_),                            \
        .nbuffers  = (_nbuffers_),                                      \
        .name      = #_name_  ,                                         \
        .model     = (&_perfmodel_)                                     \
    };

#if defined(PASTIX_STARPU_SIMULATION)

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( _name_, _nbuffers_, (starpu_cpu_func_t) 1, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( _name_, _nbuffers_, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#else

#define CODELETS_CPU(_name_, _nbuffers_ )                                  \
    CODELETS_ALL( _name_, _nbuffers_, cl_##_name_##_cpu, NULL, STARPU_CPU, 0 )

#define CODELETS_CPU_MODEL(_name_, _nbuffers_ ,_perfmodel_)                                  \
    CODELETS_ALL_MODEL( _name_, _nbuffers_, cl_##_name_##_cpu, NULL, STARPU_CPU, 0, _perfmodel_  )

#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)                       \
    CODELETS_ALL( _name_, _nbuffers_, cl_##_name_##_cpu, cl_##_name_##_gpu, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )

#define CODELETS_GPU_MODEL(_name_, _nbuffers_, _cuda_flags_, _perfmodel_)                       \
    CODELETS_ALL_MODEL(  _name_, _nbuffers_,cl_##_name_##_cpu, cl_##_name_##_gpu, STARPU_CPU | STARPU_CUDA, _cuda_flags_ , _perfmodel_ )
#endif

#if !defined(PASTIX_WITH_CUDA)
#undef CODELETS_GPU
#define CODELETS_GPU(_name_, _nbuffers_, _cuda_flags_)       \
    CODELETS_CPU( _name_, _nbuffers_ )
#endif


#endif /* _codelets_h_ */
