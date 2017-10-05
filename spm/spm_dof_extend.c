/**
 *
 * @file spm_dof_extend.c
 *
 * SParse Matrix package random multi-dofs generator.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Alban Bellot
 * @date 2015-01-01
 *
 **/
#include "common.h"
#include "spm.h"

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Generate a random multidof spm from a given spm (with dof=1).
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix used to generate the new multidof spm.
 *
 * @param[in] type
 *          Defines how to generate dofs.
 *          - 0: Generate a constant dof vector,
 *          - else: Generate a variable dof vector.
 *
 * @param[in] dof
 *          The maximum value for dofs.
 *
 ********************************************************************************
 *
 * @return the new multidof spm.
 *
 *******************************************************************************/
pastix_spm_t *
spmDofExtend( const pastix_spm_t *spm,
              const int           type,
              const int           dof )
{
    pastix_spm_t *newspm;

    /* Quick return */
    if ( dof == 1 )
        return (pastix_spm_t *)spm;

    if ( spm->dof != 1 ) {
        pastix_print_error( "Cannot extend spm including dofs already\n" );
        return (pastix_spm_t *)spm;
    }

    newspm = spmCopy( spm );

    /*
     * Generate constant dof
     */
    if (type == 0) {
        newspm->dof = dof;
    }
    else {
        pastix_int_t i, dofi, baseval;
        pastix_int_t *dofptr;

        baseval = spmFindBase( spm );

        newspm->dof  = -1;
        newspm->dofs = malloc( (spm->n+1) * sizeof(pastix_int_t) );
        dofptr = newspm->dofs;

        /*
         * Initialize the dofs array where the degree of freedom of vertex i is
         * dof[i+1] - dof[i]
         */
        *dofptr = baseval;
        for(i=0; i<spm->n; i++, dofptr++) {
            dofi = 1 + ( rand() % dof );
            dofptr[1] = dofptr[0] + dofi;
        }
    }

    spmUpdateComputedFields( newspm );

    switch (spm->flttype) {
    case PastixFloat:
        s_spmDofExtend( newspm );
        break;

    case PastixDouble:
        d_spmDofExtend( newspm );
        break;

    case PastixComplex32:
        c_spmDofExtend( newspm );
        break;

    case PastixComplex64:
        z_spmDofExtend( newspm );
        break;

    case PastixPattern:
        ;
    }

    return newspm;
}
