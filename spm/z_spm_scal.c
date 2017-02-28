/**
 * @file z_spm_scal.c
 *
 *  PaStiX spm computational routines.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
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
 * @ingroup pastix_spm_dev
 *
 * @brief Scal the spm: A = alpha * A
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *           The scaling parameter.
 *
 * @param[in,out] spm
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
