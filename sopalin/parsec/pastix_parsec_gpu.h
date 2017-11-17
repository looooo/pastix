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

#include <parsec.h>
#include <parsec/devices/device.h>

extern volatile int32_t *parsec_nbtasks_on_gpu;

/**
 *
 * TODO: Let's not forget to add some documentation in the final version
 */
static inline int
pastix_parsec_selectgpu_fct( const void *arg,
                             double      weight )
{
    (void)arg;
    (void)weight;
#if defined(PASTIX_GENERATE_MODEL)
    /* { */
    /*     static int dev_id = -1; */

    /*     return (dev_id++) % (parsec_devices_enabled()-2); */
    /* } */
    return 0;
#else
    return -2;
#endif
}

/**
 *
 * TODO: Let's not forget to add some documentation in the final version
 */
int
pastix_parsec_selectgpu_gemm1d_3p_fct( void    *arg,
                                       double   ratio,
                                       double   cpu_cost,
                                       double   gpu_cost,
                                       pastix_int_t brownum );

#endif /* _pastix_parsec_gpu_h_ */

/**
 *@}
 */
