/**
 *
 * @file pastix_parsec_gpu.h
 *
 * PaRSEC GPU functions for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @addtogroup pastix_parsec
 * @{
 *
 **/
#ifndef _pastix_parsec_gpu_h_
#define _pastix_parsec_gpu_h_

static inline int
pastix_parsec_selectgpu_fct( const void *arg,
                             double      weight )
{
    (void)arg;
    (void)weight;
#if defined(PASTIX_GENERATE_MODEL)
    return 0;
#else
    return -2;
#endif
}

#endif /* _pastix_parsec_gpu_h_ */

/**
 *@}
 */
