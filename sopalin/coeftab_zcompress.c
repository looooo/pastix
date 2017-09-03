/**
 *
 * @file coeftab_zcompress.c
 *
 * Precision dependent function to compress/uncompress the coefficients
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2017-04-28
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include <lapacke.h>
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/**
 *******************************************************************************
 *
 * @brief Compress all the cblks marked as valid for low-rank format.
 *
 * All the cblk in the top levels of the elimination tree markes as candidates
 * for compression are compressed if there is a gain to compress them. The
 * compression to low-rank format is parameterized by the input information
 * stored in the lowrank structure. On exit, all the cblks marked for
 * compression are stored through the low-rank structure, even if they are kept
 * in their full-rank form.
 *
 * @remark This routine is sequential
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem to compress.
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 * Bytes.
 *
 *******************************************************************************/
pastix_int_t
coeftab_zcompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    pastix_int_t cblknum, gain = 0;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            gain += cpucblk_zcompress( side, cblk, solvmtx->lowrank );
        }
    }
    return gain;
}

/**
 *******************************************************************************
 *
 * @brief Uncompress all column block in low-rank format into full-rank format
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem.
 *
 *******************************************************************************/
void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk  *cblk   = solvmtx->cblktab;
    pastix_int_t cblknum;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (cblk->cblktype & CBLK_COMPRESSED) {
            cpucblk_zuncompress( side, cblk );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the memory gain of the low-rank form over the full-rank form
 * for the entire matrix.
 *
 * This function returns the memory gain in bytes for the full matrix when
 * column blocks are stored in low-rank format compared to a full rank storage.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix of the problem.
 *
 *******************************************************************************
 *
 * @return The difference in favor of the low-rank storage against the full rank
 *         storage.
 *
 *******************************************************************************/
pastix_int_t
coeftab_zmemory( const SolverMatrix *solvmtx )
{
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    SolverCblk  *cblk = solvmtx->cblktab;
    pastix_int_t cblknum;
    pastix_int_t gain = 0;
    pastix_int_t original = 0;
    double       memgain, memoriginal;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        original += cblk_colnbr( cblk ) * cblk->stride;
        if (cblk->cblktype & CBLK_COMPRESSED) {
            gain += cpucblk_zmemory( side, cblk );
        }
    }

    if ( side == PastixLUCoef ) {
        original *= 2;
    }

    memgain     = gain     * pastix_size_of( PastixComplex64 );
    memoriginal = original * pastix_size_of( PastixComplex64 );
    pastix_print(0, 0,
                 OUT_LOWRANK_SUMMARY,
                 (long)gain, (long)original,
                 MEMORY_WRITE(memgain),     MEMORY_UNIT_WRITE(memgain),
                 MEMORY_WRITE(memoriginal), MEMORY_UNIT_WRITE(memoriginal));

    return gain;
}
