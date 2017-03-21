/**
 * @file z_spm_scal.c
 *
 * SParse Matrix package scaling routine.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-06-01
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_scal
 *
 * @brief Scal the spm: A = alpha * A
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[inout] spm
 *           The spm which needs to be scaled.
 *
 *******************************************************************************/
void
z_spmScal( const pastix_complex64_t alpha, pastix_spm_t *spm )
{
    pastix_int_t        nnz, i;
    pastix_complex64_t *values;

    nnz    = spm->nnz;
    values = spm->values;

    for (i=0; i<nnz; i++){
        values[i] *= alpha;
    }
}
