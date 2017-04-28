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
#include "coeftab.h"
#include "pastix_zcores.h"

/* Section: Functions */
pastix_int_t
coeftab_zcompress_one( SolverCblk *cblk,
                       pastix_lr_t lowrank )
{
    pastix_lrblock_t   *LRblocks;
    SolverBlok         *blok     = cblk[0].fblokptr;
    SolverBlok         *lblok    = cblk[1].fblokptr;
    pastix_complex64_t *lcoeftab = cblk->lcoeftab;
    pastix_complex64_t *ucoeftab = cblk->ucoeftab;
    pastix_int_t        ncols    = cblk_colnbr( cblk );
    pastix_int_t        gainL    = ncols * cblk->stride;
    pastix_int_t        gainU    = ncols * cblk->stride;
    int factoLU = (cblk->ucoeftab == NULL) ? 0 : 1;

    /* One allocation per cblk */
    LRblocks = malloc( (factoLU+1) * (lblok - blok) * sizeof(pastix_lrblock_t) );

    /* H then split */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    /**
     * Diagonal block (Not compressed)
     */
    core_zlralloc( ncols, ncols, -1, LRblocks );
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                         lcoeftab,    ncols,
                         LRblocks->u, LRblocks->rkmax );
    blok->LRblock = LRblocks; LRblocks++;

    gainL -= ncols * ncols;

    if (factoLU) {
        core_zlralloc( ncols, ncols, -1, LRblocks );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                             ucoeftab,    ncols,
                             LRblocks->u, LRblocks->rkmax );
        LRblocks++;

        gainU -= ncols * ncols;
    }

    for (blok++; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        blok->LRblock = LRblocks;

        if ( (ncols > lowrank.compress_min_width ) &&
             (nrows > lowrank.compress_min_height) )
        {
            lowrank.core_ge2lr( lowrank.tolerance, nrows, ncols,
                                lcoeftab + blok->coefind, nrows,
                                blok->LRblock );
            gainL -= (LRblocks->rk == -1) ? nrows * ncols
                : ((nrows+ncols) * LRblocks->rk);
        }
        else
        {
            core_zlralloc( nrows, ncols, -1, LRblocks );
            LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', nrows, ncols,
                                 lcoeftab + blok->coefind, nrows,
                                 LRblocks->u, LRblocks->rkmax );
        }

        LRblocks++;

        if (factoLU) {
            if ( (ncols > lowrank.compress_min_width ) &&
                 (nrows > lowrank.compress_min_height) )
            {
                lowrank.core_ge2lr( lowrank.tolerance, nrows, ncols,
                                    ucoeftab + blok->coefind, nrows,
                                    blok->LRblock+1 );
                gainU -= (LRblocks->rk == -1) ? nrows * ncols
                    : ((nrows+ncols) * LRblocks->rk);
            }
            else
            {
                core_zlralloc( nrows, ncols, -1, LRblocks );
                LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', nrows, ncols,
                                     ucoeftab + blok->coefind, nrows,
                                     LRblocks->u, LRblocks->rkmax );
            }

            LRblocks++;
        }
    }

    cblk->cblktype |= CBLK_COMPRESSED;
    /**
     * Free the dense version
     */
    free(cblk->lcoeftab); cblk->lcoeftab = NULL;
    if (cblk->ucoeftab) {
        free(cblk->ucoeftab); cblk->ucoeftab = NULL;
    }

    return gainL + gainU;
}

void
coeftab_zalloc_one( SolverCblk *cblk )
{
    pastix_lrblock_t   *LRblocks;
    SolverBlok         *blok     = cblk[0].fblokptr;
    SolverBlok         *lblok    = cblk[1].fblokptr;
    pastix_complex64_t *lcoeftab = cblk->lcoeftab;
    pastix_complex64_t *ucoeftab = cblk->ucoeftab;
    pastix_int_t        ncols    = cblk_colnbr( cblk );
    int factoLU = (cblk->ucoeftab == NULL) ? 0 : 1;

    /* One allocation per cblk */
    LRblocks = malloc( (factoLU+1) * (lblok - blok) * sizeof(pastix_lrblock_t) );

    /* H then split */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    /**
     * Diagonal block (Not compressed)
     */
    core_zlralloc( ncols, ncols, -1, LRblocks );
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                         lcoeftab,    ncols,
                         LRblocks->u, LRblocks->rkmax );
    blok->LRblock = LRblocks; LRblocks++;

    if (factoLU) {
        core_zlralloc( ncols, ncols, -1, LRblocks );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                             ucoeftab,    ncols,
                             LRblocks->u, LRblocks->rkmax );
        LRblocks++;
    }
    for (blok++; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        blok->LRblock = LRblocks;

        core_zlralloc( nrows, ncols, -1, LRblocks );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', nrows, ncols,
                             lcoeftab + blok->coefind, nrows,
                             LRblocks->u, LRblocks->rkmax );
        LRblocks++;

        if (factoLU) {

            core_zlralloc( nrows, ncols, -1, LRblocks );
            LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', nrows, ncols,
                                 ucoeftab + blok->coefind, nrows,
                                 LRblocks->u, LRblocks->rkmax );
            LRblocks++;
        }
    }

    cblk->cblktype |= CBLK_COMPRESSED;
    /**
     * Free the dense version
     */
    free(cblk->lcoeftab); cblk->lcoeftab = NULL;
    if (cblk->ucoeftab) {
        free(cblk->ucoeftab); cblk->ucoeftab = NULL;
    }
}

pastix_int_t
coeftab_zmemory_one( SolverCblk *cblk, int factoLU )
{
    SolverBlok *blok  = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_int_t gainL = 0;
    pastix_int_t gainU = 0;

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        if (blok->LRblock[0].rk >= 0){
            gainL += ((nrows * ncols) - ((nrows+ncols) * blok->LRblock[0].rkmax));
        }

        if (factoLU) {
            if (blok->LRblock[1].rk >= 0){
                gainU += ((nrows * ncols) - ((nrows+ncols) * blok->LRblock[1].rkmax));
            }
        }
    }

    return gainL + gainU;
}

void
coeftab_zuncompress_one( SolverCblk *cblk, int factoLU )
{
    SolverBlok *blok  = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_complex64_t *lcoeftab = NULL;
    pastix_complex64_t *ucoeftab = NULL;
    int ret;

    /* One allocation per cblk */
    assert( cblk->lcoeftab == NULL );
    lcoeftab = malloc( cblk->stride * ncols * sizeof(pastix_complex64_t) );

    if ( factoLU ) {
        assert( cblk->ucoeftab == NULL );
        ucoeftab = malloc( cblk->stride * ncols * sizeof(pastix_complex64_t) );
    }

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                           blok->LRblock,
                           lcoeftab + blok->coefind, nrows );
        assert( ret == 0 );
        core_zlrfree( blok->LRblock );

        if (factoLU) {
            ret = core_zlr2ge( PastixNoTrans, nrows, ncols,
                               blok->LRblock+1,
                               ucoeftab + blok->coefind, nrows );
            assert( ret == 0 );
            core_zlrfree( blok->LRblock+1 );
        }
    }

    cblk->lcoeftab = lcoeftab;
    cblk->ucoeftab = ucoeftab;

    /**
     * Free all the LRblock structures associated to the cblk
     */
    cblk->cblktype &= ~(CBLK_COMPRESSED);
    free(cblk->fblokptr->LRblock);

    (void)ret;
}


void
coeftab_zmemory( SolverMatrix *solvmtx )
{
    SolverCblk  *cblk   = solvmtx->cblktab;
    pastix_int_t cblknum;
    int          factoLU = (solvmtx->factotype == PastixFactLU) ? 1 : 0;
    pastix_int_t gain = 0;
    pastix_int_t original = 0;
    double       memgain, memoriginal;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        original += cblk_colnbr( cblk ) * cblk->stride;
        if (cblk->cblktype & CBLK_COMPRESSED) {
            gain += coeftab_zmemory_one( cblk, factoLU );
        }
    }

    if ( factoLU ) {
        original *= 2;
    }

    memgain     = gain     * pastix_size_of( PastixComplex64 );
    memoriginal = original * pastix_size_of( PastixComplex64 );
    pastix_print(0, 0,
                 OUT_LOWRANK_SUMMARY,
                 gain, original,
                 MEMORY_WRITE(memgain),     MEMORY_UNIT_WRITE(memgain),
                 MEMORY_WRITE(memoriginal), MEMORY_UNIT_WRITE(memoriginal));
}

void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk  *cblk   = solvmtx->cblktab;
    pastix_int_t cblknum;
    int          factoLU = (solvmtx->factotype == PastixFactLU) ? 1 : 0;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (cblk->cblktype & CBLK_COMPRESSED) {
            coeftab_zuncompress_one( cblk, factoLU );
        }
    }
}

void
coeftab_zcompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_int_t cblknum;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (!(cblk->cblktype & CBLK_COMPRESSED)) {
            coeftab_zcompress_one( cblk, solvmtx->lowrank );
        }
    }
}
