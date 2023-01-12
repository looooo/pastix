/**
 *
 * @file coeftab_zcinit.c
 *
 * Mixed-Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Brieuc Nicolas
 * @date 2022-10-18
 *
 * @precisions mixed zc -> ds
 *
 **/
#define _GNU_SOURCE 1
#include "common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "sopalin/coeftab.h"
#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "pastix_zccores.h"
#include "pastix_ccores.h"

/**
 *******************************************************************************
 *
 * @brief Fully initialize a single mixed-precision cblk.
 *
 * The cblk is allocated, intialized from the bcsc, and compressed if necessary.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          The solver matrix data structure.
 *
 * @param[in] bcsc
 *          The internal block CSC structure to fill-in the matrix.
 *
 * @param[in] itercblk
 *          The index of the cblk to initialize.
 *
 * @param[inout] directory
 *          The pointer to the temporary directory where to store the output
 *          files.  Used only if PASTIX_DEBUG_DUMP_COEFTAB is defined.
 *
 *******************************************************************************/
void
cpucblk_zcinit( pastix_coefside_t    side,
                const SolverMatrix  *solvmtx,
                const pastix_bcsc_t *bcsc,
                pastix_int_t         itercblk,
                const char          *directory )
{
    SolverCblk  *cblk    = solvmtx->cblktab + itercblk;
    int          ilukmax = solvmtx->lowrank.ilu_lvl;
    int rc;

    /* Do not allocate if already allocated */
    if ( !solvmtx->globalalloc ) {
        cpucblk_calloc( side, cblk );
    }

    rc = cpucblk_zcfillin( side, solvmtx, bcsc, itercblk );
    if( rc != 0 ) {
        pastix_print_error("cpucblk_zcinit: mixed-precision overflow during the matrix initialization");
        return;
    }

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /*
     * Rk: This function is not in the kernel directory to avoid the double
     * dependency with the pastix library due to pastix_fopenw()
     */
    {
        cpucblk_cdumpfile(side, cblk, itercblk, directory);
    }
#endif /* defined(PASTIX_DEBUG_DUMP_COEFTAB) */

    /* Update ILU levels if needed */
    if ( (ilukmax > 0) && (ilukmax < INT_MAX) ) {
#if defined(PASTIX_WITH_MPI)
        /* For now, we can't compute ILU(k) levels in distributed */
        if ( solvmtx->clustnbr == 1 )
#endif
        {
            do { pastix_yield(); } while( cblk->ctrbcnt > 0 );
            coeftabComputeCblkILULevels( solvmtx, cblk );
        }
    }

    /**
     * Try to compress the cblk if needs to be compressed
     */
    if ( (cblk->cblktype & CBLK_COMPRESSED) &&
         (ilukmax < INT_MAX) )
    {
        cpucblk_ccompress( solvmtx, side, ilukmax, cblk );
    }

    (void)directory;
}
