/**
 *
 * @file coeftab_zdump.c
 *
 * Precision dependent routines to dump a solver matrix into a file when debuging.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2017-04-28
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin/coeftab_z.h"

/**
 *******************************************************************************
 *
 * @brief Dump a single column block into a FILE in a human readale format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          The column block to dump into the file.
 *
 * @param[in] array
 *          The generic pointer to lcoeftab or ucoeftab from the cblk, depending
 *          on wheter L or U is printed.
 *
 * @param[inout] stream
 *          The FILE structure opened in write mode.
 *
 *******************************************************************************/
void
coeftab_zdumpcblk( const SolverCblk *cblk,
                   const void       *array,
                   FILE             *stream )
{
    const pastix_complex64_t *coeftab = (const pastix_complex64_t*)array;
    SolverBlok *blok;
    pastix_int_t itercol;
    pastix_int_t iterrow;
    pastix_int_t coefindx;

    /* We don't know how to dump the compressed block for now */
    if ( cblk->cblktype & CBLK_COMPRESSED )
        return;

    for (itercol  = cblk->fcolnum;
         itercol <= cblk->lcolnum;
         itercol++)
    {
        /* Diagonal Block */
        blok     = cblk->fblokptr;
        coefindx = blok->coefind;
        if (cblk->cblktype & CBLK_LAYOUT_2D) {
            coefindx += (itercol - cblk->fcolnum) * blok_rownbr( blok );
        }
        else {
            coefindx += (itercol - cblk->fcolnum) * cblk->stride;
        }

        for (iterrow  = blok->frownum;
             iterrow <= blok->lrownum;
             iterrow++, coefindx++)
        {
            if ((cabs( coeftab[coefindx] ) > 0.) &&
                (itercol <= iterrow))
            {
#if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf(stream, "%ld %ld (%13e,%13e)\n",
                        (long)itercol, (long)iterrow,
                        creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                fprintf(stream, "%ld %ld %13e\n",
                        (long)itercol, (long)iterrow,
                        coeftab[coefindx]);
#endif
            }
        }

        /* Off diagonal blocks */
        blok++;
        while( blok < (cblk+1)->fblokptr )
        {
            coefindx  = blok->coefind;
            if (cblk->cblktype & CBLK_LAYOUT_2D) {
                coefindx += (itercol - cblk->fcolnum) * blok_rownbr( blok );
            }
            else {
                coefindx += (itercol - cblk->fcolnum) * cblk->stride;
            }

            for (iterrow  = blok->frownum;
                 iterrow <= blok->lrownum;
                 iterrow++, coefindx++)
            {
                if (cabs( coeftab[coefindx]) > 0.)
                {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e)\n",
                            (long)itercol, (long)iterrow,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e\n",
                            (long)itercol, (long)iterrow,
                            coeftab[coefindx]);
#endif
                }
            }
            blok++;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Dump the sovler matrix coefficients into a file in human readable
 * format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix to print.
 *
 * @param[in] filename
 *          The filename where to store the output matrix.
 *
 *******************************************************************************/
void
coeftab_zdump( const SolverMatrix *solvmtx,
               const char   *filename )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    FILE *stream = fopen( filename, "w" );

    /*
     * TODO: there is a problem right here for now, because there are no
     * distinctions between L and U coeffcients in the final file
     */

    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        coeftab_zdumpcblk( cblk, cblk->lcoeftab, stream );
        if ( NULL != cblk->ucoeftab )
            coeftab_zdumpcblk( cblk, cblk->ucoeftab, stream );
    }

    fclose( stream );
}
