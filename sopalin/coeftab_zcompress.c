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
    int ret;

    /* One allocation per cblk */
    LRblocks = malloc( (factoLU+1) * (lblok - blok) * sizeof(pastix_lrblock_t) );

    /**
     * Diagonal block (Not compressed)
     */
    core_zlralloc( ncols, ncols, -1, LRblocks );
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                         lcoeftab,    cblk->stride,
                         LRblocks->u, LRblocks->rkmax );
    blok->LRblock = LRblocks; LRblocks++;

    gainL -= ncols * ncols;

    if (factoLU) {
        core_zlralloc( ncols, ncols, -1, LRblocks );
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                             ucoeftab,    cblk->stride,
                             LRblocks->u, LRblocks->rkmax );
        LRblocks++;

        gainU -= ncols * ncols;
    }

    for (blok++; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        blok->LRblock = LRblocks;
        ret = core_zge2lr( tol, nrows, ncols,
                           lcoeftab + blok->coefind, cblk->stride,
                           blok->LRblock );
        assert( ret == 0 );

        gainL -= (LRblocks->rk == -1) ? nrows * ncols
            : ((nrows+ncols) * LRblocks->rk);
        LRblocks++;

        if (factoLU) {
            ret = core_zge2lr( tol, nrows, ncols,
                               ucoeftab + blok->coefind, cblk->stride,
                               blok->LRblock+1 );
            assert( ret == 0 );

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

pastix_int_t
coeftab_zuncompress_one( SolverCblk *cblk, int factoLU )
{
    SolverBlok *blok  = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_int_t        gainL    = 0;
    pastix_int_t        gainU    = 0;
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

        if (blok->LRblock[0].rk >= 0){
            gainL += ((nrows * ncols) - ((nrows+ncols) * blok->LRblock[0].rkmax));
        }
        ret = core_zlr2ge( nrows, ncols,
                           blok->LRblock,
                           lcoeftab + blok->coefind, cblk->stride );
        assert( ret == 0 );
        core_zlrfree( blok->LRblock );

        if (factoLU) {
            if (blok->LRblock[1].rk >= 0){
                gainU += ((nrows * ncols) - ((nrows+ncols) * blok->LRblock[1].rkmax));
            }
            ret = core_zlr2ge( nrows, ncols,
                               blok->LRblock+1,
                               ucoeftab + blok->coefind, cblk->stride );
            assert( ret == 0 );
            core_zlrfree( blok->LRblock+1 );
        }
    }

    cblk->lcoeftab = lcoeftab;
    cblk->ucoeftab = ucoeftab;
    cblk->cblktype |= CBLK_DENSE;

    /**
     * Free all the LRblock structures associated to the cblk
     */
    free(cblk->fblokptr->LRblock);
    return gainL + gainU;
}


void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_int_t cblknum;
    char  *tolerance = getenv("TOLERANCE");
    double tol = atof(tolerance);
    pastix_int_t gain = 0;
    pastix_int_t original = 0;
    double memgain, memoriginal;
    int factoLU = 1;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        original += cblk_colnbr( cblk ) * cblk->stride;
        if (!(cblk->cblktype & CBLK_DENSE)) {
            gain += coeftab_zuncompress_one( cblk, 1 );
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
    int factoLU = 1;

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

