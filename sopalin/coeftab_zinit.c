/**
 *
 * @file coeftab_zinit.c
 *
 * Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include "sopalin/coeftab_z.h"

void
coeftab_zffbcsc( const SolverMatrix  *solvmtx,
                 const pastix_bcsc_t *bcsc,
                 pastix_int_t         itercblk )
{
    const bcsc_format_t *csccblk = bcsc->cscftab + itercblk;
    SolverCblk *solvcblk = solvmtx->cblktab + itercblk;
    SolverBlok *solvblok;
    SolverBlok *solvblok2 = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab = solvcblk->lcoeftab;
    pastix_complex64_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues = bcsc->Lvalues;
    pastix_complex64_t *Uvalues = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx;

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];
        solvblok = solvcblk->fblokptr;

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < solvblok2) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                }

                if ( solvblok < solvblok2 )
                {
                    coefindx  = solvblok->coefind;
                    coefindx += rownum - solvblok->frownum;
                    if (solvcblk->cblktype & CBLK_LAYOUT_2D) {
                        coefindx += itercoltab * blok_rownbr( solvblok );
                    }
                    else {
                        coefindx += itercoltab * solvcblk->stride;
                    }

                    lcoeftab[coefindx] = Lvalues[iterval];

                    if ( (ucoeftab != NULL) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
                }
                else {
                    /* printf("ILU: csc2solv drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n", */
                    /*        (long)datacode->cblktab[itercblk].fcolnum+ */
                    /*        (long)itercoltab,(long)itercoltab, */
                    /*        (long)CSC_ROW(cscmtx,iterval),(long)iterval, */
                    /*        (long)itercblk, */
                    /*        (long)datacode->cblktab[itercblk].fcolnum, */
                    /*        (long)datacode->cblktab[itercblk].lcolnum); */
                }
            }
        }
    }
}

/*
 * Function: z_CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void
coeftab_zcblkfill( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t itercblk,
                   int factoLU )
{
    SolverCblk *cblk     = solvmtx->cblktab + itercblk;
    pastix_int_t coefnbr = cblk->stride * cblk_colnbr( cblk );

    /* If not NULL, allocated to store the shur complement for exemple */
    assert( cblk->lcoeftab == NULL );

    MALLOC_INTERN( cblk->lcoeftab, coefnbr, pastix_complex64_t );
    memset( cblk->lcoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );

    if ( factoLU ) {
        /* Extra diagonal block for low-rank updates */
        MALLOC_INTERN( cblk->ucoeftab, coefnbr, pastix_complex64_t );
        memset( cblk->ucoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );
    }
    else {
        cblk->ucoeftab = NULL;
    }

    coeftab_zffbcsc( solvmtx, bcsc, itercblk );
}

/*
 * Function: z_CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void
coeftab_zcblkinit( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t         itercblk,
                   int                  factoLU,
                   char               **directory )
{
    SolverCblk *cblk           = solvmtx->cblktab + itercblk;
    pastix_int_t compress_when = solvmtx->lowrank.compress_when;

    coeftab_zcblkfill( solvmtx, bcsc, itercblk, factoLU );

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    {
        FILE *f = NULL;
        char *filename;
        int rc;

        rc = asprintf( &filename, "Lcblk%05ld_init.txt", itercblk );
        f  = pastix_fopenw( directory, filename, "w" );
        if ( f != NULL ) {
            coeftab_zcblkdump( cblk, PastixLower, f );
            fclose( f );
        }
        free( filename );

        if ( cblk->ucoeftab ) {
            rc = asprintf( &filename, "Ucblk%05ld_init.txt", itercblk );
            f  = pastix_fopenw( directory, filename, "w" );
            if ( f != NULL ) {
                coeftab_zcblkdump( cblk, PastixUpper, f );
                fclose( f );
            }
            free( filename );
        }
        (void)rc;
    }
#endif /* defined(PASTIX_DEBUG_DUMP_COEFTAB) */

    /**
     * Try to compress the cblk if needs to be compressed
     * TODO: change the criteria based on the level in the tree
     */
    if ( cblk->cblktype & CBLK_LAYOUT_2D )
    {
        if ( compress_when == PastixCompressWhenBegin ){
            coeftab_zcompress_one( cblk, solvmtx->lowrank );
        }
        else if ( compress_when == PastixCompressWhenEnd ){
            coeftab_zalloc_one( cblk );
        }
    }
}
