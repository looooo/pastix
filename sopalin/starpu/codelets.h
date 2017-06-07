/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
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

#if defined(PASTIX_STARPU_SIMULATION)

#define CODELETS_ALL( name, _nbuffers, _original_location_, cuda_flags ) \
    struct starpu_codelet cl_##name = {                                 \
        .where     = STARPU_CPU,                                        \
        .cpu_func  = (starpu_cpu_func_t) 1,                             \
        .cuda_func = NULL,                                              \
        .nbuffers  = (_nbuffers),                                       \
        .name      = "cl_" #name                                        \
    };

#else

#define CODELETS_ALL( name, _nbuffers, _original_location_, cuda_flags ) \
    struct starpu_codelet cl_##name = {                                 \
        .where     = (_original_location_),                             \
        .cpu_func  = (cl_##name##_cpu),                                 \
        CODELET_CUDA_FLAGS(cuda_flags)                                  \
        .cuda_func = (cl_##name##_gpu),                                 \
        .nbuffers  = (_nbuffers),                                       \
        .name      = "cl_" #name                                        \
    };

#endif

#define CODELETS_CPU(name, _nbuffers )                  \
    CODELETS_ALL( name, _nbuffers, STARPU_CPU, 0 )

#if defined(PASTIX_WITH_CUDA)
#define CODELETS_GPU(name, _nbuffers, cuda_flags)                       \
    CODELETS_ALL( name, _nbuffers, STARPU_CPU  | STARPU_CUDA, cuda_flags )
#else
#define CODELETS_GPU(name, _nbuffers, cuda_flags)       \
    CODELETS_CPU( name, _nbuffers )
#endif

#endif /* _CODELETS_H_ */
