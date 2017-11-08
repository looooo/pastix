/**
 *
 * @file models.h
 *
 * PaStiX performance models routines
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2017-11-09
 *
 * @addtogroup pastix_models
 * @{
 *    This module contains all the function to load performance models that will
 *    be used afterwards by simulation, and runtime schedulers.
 *
 **/
#ifndef _models_h_
#define _models_h_

#include "kernels/kernels_trace.h"

typedef struct pastix_model_s {
    char *name;
    double coefficients[4][PastixKernelLvl1Nbr][8];
} pastix_model_t;

static inline double
modelsGetCost1Param( const double *coefs, pastix_int_t N )
{
    /* a3 * N^3 + a2 * N^2 + a1 * N + a0 */
    return ((coefs[3] * N + coefs[2]) * N + coefs[1]) * N + coefs[0];
}

static inline double
modelsGetCost2Param( const double *coefs, pastix_int_t M, pastix_int_t N )
{
    /* a5 * M*N^2 + a4 * M*N + a3 * N^2 + a2 * M + a1 * N + a0 */
    return ((coefs[5] * (double)M + coefs[3]) * (double)N + coefs[4] * (double)M + coefs[1]) * (double)N + coefs[2] * (double)M + coefs[0];
}

static inline double
modelsGetCost3Param( const double *coefs, pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* a7 * M * N * K + a6 * M * K + a5 * K * N + a4 * M * N + a3 * M + a2 * N + a1 * K + a0 */
    return (coefs[7] * (double)M * (double)N * (double)K +
            coefs[6] * (double)M * (double)K +
            coefs[5] * (double)K * (double)N +
            coefs[4] * (double)M * (double)N +
            coefs[3] * (double)M +
            coefs[2] * (double)N +
            coefs[1] * (double)K +
            coefs[0]);
}

void modelsLoad( pastix_data_t *pastix_data );

#endif /* _models_h_ */

/**
 * @}
 */
