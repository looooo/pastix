/**
 *
 * @file coeftab_zinit.c
 *
 * Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @date 2023-07-21
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "sopalin/coeftab.h"
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
 * @param[in] side
 *          Define which side of the cblk must be printed.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *
 * @param[in] cblk
 *          The column block to dump into the file.
 *
 * @param[inout] stream
 *          The FILE structure opened in write mode.
 *
 *******************************************************************************/
void
cpucblk_zdump( pastix_coefside_t side,
               const SolverCblk *cblk,
               FILE             *stream )
{
    const pastix_complex64_t *coeftab = side == PastixUCoef ? cblk->ucoeftab : cblk->lcoeftab;
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
                if ( side == PastixUCoef ) {
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
                    if ( side == PastixUCoef ) {
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
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          TODO
 *
 * @param[in] cblk
 *          TODO
 *
 *******************************************************************************/
void coeftab_zcblkComputePreselect( const SolverMatrix *solvmtx, SolverCblk *cblk )
{
    /* If there are off diagonal supernodes in the column */
    SolverBlok *blok = cblk->fblokptr + 1; /* this diagonal block */
    SolverBlok *lblk = cblk[1].fblokptr;   /* the next diagonal block */

    for (; blok<lblk; blok++) {
        SolverCblk *fcblk = solvmtx->cblktab + blok->fcblknm;
        int is_preselected = blok_is_preselected( cblk, blok, fcblk );

        if ( is_preselected ) {
            blok->iluklvl = -1;
        }
        else {
            blok->iluklvl = INT_MAX;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Fully initialize a single cblk.
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
cpucblk_zinit( pastix_coefside_t    side,
               const SolverMatrix  *solvmtx,
               const pastix_bcsc_t *bcsc,
               pastix_int_t         itercblk,
               const char          *directory )
{
    SolverCblk  *cblk    = solvmtx->cblktab + itercblk;
    int          ilukmax = solvmtx->lowrank.ilu_lvl;

    /* Do not allocate if already allocated */
    if ( !solvmtx->globalalloc ) {
        cpucblk_zalloc( side, cblk );
    }

    cpucblk_zfillin( side, solvmtx, bcsc, itercblk );

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /*
     * Rk: This function is not in the kernel directory to avoid the double
     * dependency with the pastix library due to pastix_fopenw()
     */
    {
        cpucblk_zdumpfile(side, cblk, itercblk, directory);
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
        cpucblk_zcompress( solvmtx, side, ilukmax, cblk );
    }

    (void)directory;
}
