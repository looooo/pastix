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
 * @ingroup pastix_spm_internal
 *
 * z_spmScal - Scal the matrix with ||A||_2
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *           The spm which needs to be scaled.
 *
 *******************************************************************************/
void
z_spmScal( pastix_spm_t *spm )
{
    double              norm;
    pastix_int_t        nnz, i;
    pastix_complex64_t *values;

    nnz    = spm->nnz;
    values = spm->values;
    norm   = z_spmNorm( PastixFrobeniusNorm, spm );

    for (i=0; i<nnz; i++){
        values[i] /= norm;
    }
}
