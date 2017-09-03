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
#include "pastix_zcores.h"

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
 * @param[in] uplo
 *          Specify if upper or lower part is printed.
 *
 * @param[inout] stream
 *          The FILE structure opened in write mode.
 *
 *******************************************************************************/
void
cpucblk_zdump( const SolverCblk *cblk,
               pastix_uplo_t     uplo,
               FILE             *stream )
{
    const pastix_complex64_t *coeftab = uplo == PastixUpper ? cblk->ucoeftab : cblk->lcoeftab;
    SolverBlok  *blok;
    pastix_int_t itercol;
    pastix_int_t iterrow;
    pastix_int_t coefindx;

    /* We don't know how to dump the compressed block for now */
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        fprintf(stderr, "coeftab_zcblkdump: Can't dump a compressed cblk\n");
        return;
    }

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
                if ( uplo == PastixUpper ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e) [U]\n",
                            (long)itercol, (long)iterrow,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e [U]\n",
                            (long)itercol, (long)iterrow,
                            coeftab[coefindx]);
#endif
                }
                else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e) [L]\n",
                            (long)iterrow, (long)itercol,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e [L]\n",
                            (long)iterrow, (long)itercol,
                            coeftab[coefindx]);
#endif
                }
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
                    if ( uplo == PastixUpper ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(stream, "%ld %ld (%13e,%13e) [U]\n",
                                (long)itercol, (long)iterrow,
                                creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                        fprintf(stream, "%ld %ld %13e [U]\n",
                                (long)itercol, (long)iterrow,
                                coeftab[coefindx]);
#endif
                    }
                    else {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf(stream, "%ld %ld (%13e,%13e) [L]\n",
                                (long)iterrow, (long)itercol,
                                creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                        fprintf(stream, "%ld %ld %13e [L]\n",
                                (long)iterrow, (long)itercol,
                                coeftab[coefindx]);
#endif
                    }
                }
            }
            blok++;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Dump the solver matrix coefficients into a file in human readable
 * format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data instance to access the unique directory id in which
 *          output the files.
 *
 * @param[in] solvmtx
 *          The solver matrix to print.
 *
 * @param[in] filename
 *          The filename where to store the output matrix.
 *
 *******************************************************************************/
void
coeftab_zdump( pastix_data_t      *pastix_data,
               const SolverMatrix *solvmtx,
               const char         *filename )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    FILE *stream = NULL;

    stream = pastix_fopenw( &(pastix_data->dirtemp), filename, "w" );
    if ( stream == NULL ){
        return;
    }

    /*
     * TODO: there is a problem right here for now, because there are no
     * distinctions between L and U coeffcients in the final file
     */
    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        cpucblk_zdump( cblk, PastixLower, stream );
        if ( NULL != cblk->ucoeftab )
            cpucblk_zdump( cblk, PastixUpper, stream );
    }

    fclose( stream );
}
