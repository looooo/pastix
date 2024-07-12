/**
 *
 * @file codelets.h
 *
 * @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2024-07-05
 *
 * Macro helper to define the static codelet structures used to submit the tasks.
 *
 **/
#ifndef _codelets_h_
#define _codelets_h_

#include "common.h"

/* Define a macro to set the cuda flags if defined in StarPU */
#if defined(STARPU_CUDA_ASYNC)
#define CODELET_CUDA_FLAGS( _flags_ ) .cuda_flags = {(_flags_)},
#else
#define CODELET_CUDA_FLAGS( _flags_ )
#endif

/* Skeleton macro to define the codelet structure */
#define CODELETS_SKELETON( _unit_, _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ ) \
    struct starpu_codelet cl_##_name_##_##_unit_ = {                    \
        .where     = (_original_location_),                             \
        .cpu_funcs[0] = (_cpu_func_name_),                              \
        CODELET_CUDA_FLAGS(_cuda_flags_)                                \
        .cuda_funcs[0] = (_cuda_func_name_),                            \
        .nbuffers  = (_nbuffers_),                                      \
        .name      = #_name_ ,                                          \
        .model     = &(starpu_##_name_##_model)                         \
    }


/* Add a first layer above the skeleton in case we are in a simulated run with StarPU-Simgrid */
#if defined(PASTIX_STARPU_SIMULATION)
#define CODELETS_ANY( _unit_, _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ ) \
    CODELETS_SKELETON( _unit_, _name_, _nbuffers_, (starpu_cpu_func_t) 1, (starpu_cuda_func_t) 1, _original_location_, _cuda_flags_ )
#else
#define CODELETS_ANY( _unit_, _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ ) \
    CODELETS_SKELETON( _unit_, _name_, _nbuffers_, _cpu_func_name_, _cuda_func_name_, _original_location_, _cuda_flags_ )
#endif

/* Define the CPU only codelets */
#define CODELETS_CPU( _name_, _nbuffers_ )                              \
    CODELETS_ANY( cpu, _name_, _nbuffers_, fct_##_name_##_cpu, NULL, STARPU_CPU, 0 )

/* Define the CPU/GPU codelets */
#if defined (PASTIX_WITH_CUDA)
#define CODELETS_GPU( _name_, _nbuffers_, _cuda_flags_ )                \
    CODELETS_ANY( any, _name_, _nbuffers_, fct_##_name_##_cpu, fct_##_name_##_gpu, STARPU_CPU | STARPU_CUDA, _cuda_flags_ )
#else
#define CODELETS_GPU( _name_, _nbuffers_, _cuda_flags_ )                \
    CODELETS_ANY( any, _name_, _nbuffers_, fct_##_name_##_cpu, NULL, STARPU_CPU, 0 )
#endif

#endif /* _codelets_h_ */
