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

#endif /* _pastix_parsec_gpu_h_ */

/**
 *@}
 */
