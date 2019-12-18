/**
 *
 * @file solver_check.c
 *
 * PaStiX check function fo rthe solver structure.
 *
 * @copyright 2004-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include <stdio.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "simu.h"
#include "pastix_zcores.h"
#include "pastix_ccores.h"
#include "pastix_dcores.h"
#include "pastix_scores.h"

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Checks the consistency of the given solver matrix structure.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure to check.
 *
 *******************************************************************************
 *
 * @retval 0 if the structure is correct
 * @retval 1 if incorrect
 *
 *******************************************************************************/
int
solverCheck( const SolverMatrix *solvmtx )
{
    int          i, j;
    SolverBlok  *blok, *fblok;
    SolverCblk  *cblk, *fcblk;
    pastix_int_t bloknum;

    assert( solvmtx->cblknbr <= solvmtx->gcblknbr );

    cblk = solvmtx->cblktab;

    for ( i = 0; i < solvmtx->cblknbr; i++, cblk++ ) {
        /* Make sure the lock is initialized correctly */
        assert( cblk->lock == PASTIX_ATOMIC_UNLOCKED );

        /*
         * Check that:
         *    - the brownum is lower or equal to the next one
         *    - the brown2d is included between the two brow values
         */
        assert( cblk[1].brownum >= cblk[0].brownum );
        assert( ( cblk[0].brownum <= cblk[0].brown2d ) && ( cblk[0].brown2d <= cblk[1].brownum ) );

        /* Check dimensions */
        assert( cblk->fcolnum <= cblk->lcolnum );
        assert( cblk->stride >= cblk_colnbr( cblk ) );

        fblok = cblk->fblokptr;

        assert( cblk[0].lcolnum < cblk[1].fcolnum );
        assert( fblok->browind == -1 );
        assert( fblok->lcblknm == i );

        if ( cblk->cblktype & CBLK_FANIN ) {
            /* Dimensions are similar to any other cblk but storeed remotely */
            assert( cblk->ownerid != solvmtx->clustnum );

            assert( fblok->fcblknm == -1 );
        }
        else {
            assert( (cblk[1].cblktype & CBLK_FANIN) ||
                    (!(cblk[1].cblktype & CBLK_FANIN) && (cblk[0].lcolidx < cblk[1].lcolidx)) );
            assert( (solvmtx->gcbl2loc == NULL) ||
                    (cblk->ownerid == solvmtx->clustnum) );

            assert( fblok->lcblknm == fblok->fcblknm );
        }

        /* Check bloks in the current cblk */
        blok = fblok + 1;
        for ( ; blok < cblk[1].fblokptr; blok++ ) {
            assert( blok->lcblknm == i );
            assert( blok->frownum <= blok->lrownum );

            /* Next block in the same cblk must be after */
            assert( (blok[1].lcblknm != blok[0].lcblknm) ||
                    ((blok[1].lcblknm == blok[0].lcblknm) && (blok[0].lrownum < blok[1].frownum)) );

            if ( cblk->cblktype & CBLK_FANIN ) {
                assert( blok->fcblknm == -1 );
            }
            else {
                fcblk = solvmtx->cblktab + blok->fcblknm;
                fblok = fcblk->fblokptr;

                assert( blok->lcblknm < blok->fcblknm );
                assert( blok->frownum >= fblok->frownum );
                assert( blok->lrownum <= fblok->lrownum );
            }
        }

        /* Check previous bloks in row */
        for ( j = cblk[0].brownum; j < cblk[1].brownum; j++ ) {
            bloknum = solvmtx->browtab[j];
            blok    = solvmtx->bloktab + bloknum;
            fcblk   = solvmtx->cblktab + blok->lcblknm;

            assert( blok->browind == j );
            assert( blok->fcblknm == i );

            /* Bloks in Fanin sould never appear in the browtab */
            assert( !(fcblk->cblktype & CBLK_FANIN) );
            assert( fcblk->ownerid == solvmtx->clustnum );
        }
    }

    return 0;
}

