/**
 *
 * @file coeftab_zcompress.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
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
                       double      tol )
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

    /**
     * Diagonal block (Not compressed)
     */
    LRblocks->rk    = -1;
    LRblocks->rkmax = -1;
    LRblocks->u     = malloc( ncols * ncols * sizeof(pastix_complex64_t) );
    LRblocks->v     = NULL;
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                         lcoeftab, cblk->stride, LRblocks->u, ncols );
    blok->coefL_u_LR = LRblocks->u;
    blok->coefL_v_LR = LRblocks->v;
    blok->rankL      = LRblocks->rk;
    blok->LRblock = LRblocks; LRblocks++;

    gainL -= ncols * ncols;

    if (factoLU) {
        LRblocks->rk    = -1;
        LRblocks->rkmax = -1;
        LRblocks->u     = malloc( ncols * ncols * sizeof(pastix_complex64_t) );
        LRblocks->v     = NULL;
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                             ucoeftab, cblk->stride, LRblocks->u, ncols );
        blok->coefU_u_LR = LRblocks->u;
        blok->coefU_v_LR = LRblocks->v;
        blok->rankU      = LRblocks->rk;
        LRblocks++;

        gainU -= ncols * ncols;
    }

    for (blok++; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        blok->LRblock = LRblocks;
        core_zge2lr( tol, nrows, ncols,
                     lcoeftab + blok->coefind, cblk->stride,
                     blok->LRblock );
        blok->coefL_u_LR = LRblocks->u;
        blok->coefL_v_LR = LRblocks->v;
        blok->rankL      = LRblocks->rk;

        gainL -= (LRblocks->rk == -1) ? nrows * ncols
            : ((nrows+ncols) * LRblocks->rk);
        LRblocks++;

        if (factoLU) {
            core_zge2lr( tol, nrows, ncols,
                         ucoeftab + blok->coefind, cblk->stride,
                         blok->LRblock+1 );
            blok->coefU_u_LR = LRblocks->u;
            blok->coefU_v_LR = LRblocks->v;
            blok->rankU      = LRblocks->rk;
            gainU -= (LRblocks->rk == -1) ? nrows * ncols
                : ((nrows+ncols) * LRblocks->rk);
            LRblocks++;
        }
    }

    cblk->cblktype &= ~(CBLK_DENSE);
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
coeftab_zuncompress_one( SolverCblk *cblk, int factoLU )
{
    SolverBlok *blok  = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_complex64_t *lcoeftab = NULL;
    pastix_complex64_t *ucoeftab = NULL;

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

        core_zlr2ge( nrows, ncols,
                     blok->LRblock,
                     lcoeftab + blok->coefind, cblk->stride );

        free( blok->LRblock[0].u ); blok->LRblock[0].u = NULL;
        free( blok->LRblock[0].v ); blok->LRblock[0].v = NULL;
        if (factoLU) {
            core_zlr2ge( nrows, ncols,
                         blok->LRblock+1,
                         ucoeftab + blok->coefind, cblk->stride );
            free( blok->LRblock[1].u ); blok->LRblock[1].u = NULL;
            free( blok->LRblock[1].v ); blok->LRblock[1].v = NULL;
       }
    }

    cblk->lcoeftab = lcoeftab;
    cblk->ucoeftab = ucoeftab;
    cblk->cblktype |= CBLK_DENSE;

    /**
     * Free all the LRblock structures associated to the cblk
     */
    free(cblk->fblokptr->LRblock);
}


void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_int_t cblknum;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (!(cblk->cblktype & CBLK_DENSE)) {
            coeftab_zuncompress_one( cblk, 1 );
        }
    }
}

void
coeftab_zcompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_int_t cblknum;
    char  *tolerance = getenv("TOLERANCE");
    double tol = atof(tolerance);
    pastix_int_t gain = 0;
    pastix_int_t original = 0;
    double memgain, memoriginal;
    int factoLU = 0;

    if ( solvmtx->cblktab[0].ucoeftab ) {
        factoLU = 1;
    }

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        original += cblk_colnbr( cblk ) * cblk->stride;
        if (cblk->cblktype & CBLK_DENSE) {
            gain += coeftab_zcompress_one( cblk, tol );
        }
    }

    if ( factoLU ) {
        original *= 2;
    }

    memgain     = gain     * pastix_size_of( PastixComplex64 );
    memoriginal = original * pastix_size_of( PastixComplex64 );
    fprintf( stdout,
             "  Compression : Tolerance         %e\n"
             "                Elements removed  %ld / %ld\n"
             "                Memory saved      %.3g %s / %.3g %s\n",
             tol, gain, original,
             MEMORY_WRITE(memgain),     MEMORY_UNIT_WRITE(memgain),
             MEMORY_WRITE(memoriginal), MEMORY_UNIT_WRITE(memoriginal));
}

