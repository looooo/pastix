#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "cost.h"
#include "symbol.h"
#include "elimin.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "solver.h"
#include "simu.h"

#include "extracblk.h"

static inline
int pastix_blend_with_constant_split() {
    return pastix_env_is_set_to("PASTIX_BLEND_SPLIT", "CONSTANT");
}

static inline
int pastix_blend_with_smallest_upper_split() {
    return pastix_env_is_set_to("PASTIX_BLEND_SPLIT", "UPPER");
}

static inline
int pastix_blend_split_percent() {
    return
      pastix_getenv_get_value_int( "PASTIX_BLEND_SPLIT_AUTORIZED_PERCENTAGE",
                                    10);
}


static inline pastix_int_t
computeNbSplit( const BlendCtrl *ctrl,
                pastix_int_t     candnbr,
                pastix_int_t     width )
{
    pastix_int_t blas_min_col;
    pastix_int_t blas_max_col;
    pastix_int_t step, nseq;

    /* Compute minimun broadness for splitting this cblk */
    if(candnbr > ctrl->ratiolimit)
    {
        blas_min_col = ctrl->blblokmin;
        blas_max_col = ctrl->blblokmax;
    }
    else
    {
        blas_min_col = ctrl->blcolmin;
        blas_max_col = ctrl->blcolmax;
    }

    if(candnbr == 1)
    {
        /*** Need to split big supernode because
         the diagonal block factorization is written
         in BLAS1 (due to the pivoting in LDLt and LU) ***/
        /*
         * if the column block size is small enough there is no need to
         * split it.
         */
        if( width <= blas_max_col) {
            return 1;
        }

        nseq = pastix_iceil( width, blas_max_col );
    }
    else
    {
        pastix_int_t abs = ctrl->abs;
        if(candnbr > ctrl->ratiolimit)
        {
            abs *= 2; /* Increase abs for 2D */
        }

        /* If option adaptative block size is set then compute the size of a column block */
        if(abs > 0)
        {
            step = pastix_iceil( width, (abs * candnbr) );

            step = pastix_imax(step, blas_min_col);
            step = pastix_imin(step, blas_max_col);

            /* Ceil */
            nseq = pastix_iceil( width, step );
        }
        else
        {
            nseq = pastix_iceil( width, blas_max_col );
        }
    }

    /* Make sure cblk are at least blas_min_col wide */
    if ( nseq > 1 && (width / nseq) < blas_min_col ) {
        nseq--;
    }

    return nseq;
}

static inline pastix_int_t *
computeNbBlocksPerLine( const SymbolMatrix *symbmtx,
                        pastix_int_t frowsplit )
{
    SymbolBlok   *curblok;
    pastix_int_t *nblocksperline;
    pastix_int_t  bloknum, line;
    pastix_int_t  size = symbmtx->nodenbr - frowsplit + 1;

    /*
     * Allocate the temporary buffer nblocksperline, nbblocksperline stores the
     * number of blocks that will be splitted if with split between line i and
     * i+1.
     */
    MALLOC_INTERN( nblocksperline, size, pastix_int_t );
    memset( nblocksperline, 0, size * sizeof(pastix_int_t) );

    curblok = symbmtx->bloktab;
    for(bloknum=0; bloknum<symbmtx->bloknbr; bloknum++, curblok++ )
    {
        if ( curblok->lrownum < frowsplit )
            continue;

        /*
         * For each couple of rows (i,i+1) in the block, we increment
         * the number blocks in regard of the row i
         */
        for(line = pastix_imax( curblok->frownum, frowsplit);
            line < curblok->lrownum; line++ )
        {
            nblocksperline[ line-frowsplit ]++;
        }
    }
    assert( nblocksperline[ size-1 ] == 0 );

    return nblocksperline;
}

static inline pastix_int_t
computeSmallestSplit( pastix_int_t *nblocksperline,
                      pastix_int_t step,
                      pastix_int_t max,
                      pastix_int_t authorized_percent)
{
    pastix_int_t limit = pastix_iceil( step*authorized_percent, 100 );
    pastix_int_t i, lcolnum, nbsplit;
    pastix_int_t lmin, lmax, lavg;

    if (step >= max)
        return max-1;
    assert( step > 1 );

    lavg = step - 1;
    lmin = pastix_imax( lavg - limit - 1,  0   );
    lmax = pastix_imin( lavg + limit + 1,  max );

    lcolnum = lavg;
    nbsplit = nblocksperline[ lcolnum ];

    /* Search for the minimal split */
    for(i=lavg+1; i<lmax; i++ )
    {
        if ( nblocksperline[ i ] < nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }
    for(i=lavg-1; i>lmin; i-- )
    {
        if ( nblocksperline[ i ] < nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }

    return lcolnum;
}

static inline pastix_int_t
computeSmallestSplit_max( pastix_int_t *nblocksperline,
                          pastix_int_t step,
                          pastix_int_t max,
                          pastix_int_t authorized_percent )
{
    pastix_int_t limit = pastix_iceil( step*authorized_percent, 100 );
    pastix_int_t i, lcolnum, nbsplit;
    pastix_int_t lmin, lmax, lavg;

    if (step >= max)
        return max-1;
    assert( step > 1 );

    lavg = step - 1;
    lmin = pastix_imax( lavg - limit,  1   );
    lmax = pastix_imin( lavg + limit + 1,  max );

    lcolnum = lmin;
    nbsplit = nblocksperline[ lcolnum ];

    /* Search for the minimal split */
    for(i=lmin; i<lmax; i++ )
    {
        if ( nblocksperline[ i ] <= nbsplit )
        {
            lcolnum = i;
            nbsplit = nblocksperline[ i ];
        }
    }

    return lcolnum;
}

static inline pastix_int_t
computeConstantSplit( pastix_int_t *nblocksperline,
                      pastix_int_t step,
                      pastix_int_t max,
                      pastix_int_t authorized_percent )
{
    (void)nblocksperline;
    (void)authorized_percent;
    if (step >= max)
        return max-1;
    assert( step > 1 );
    return step-1;
}

/*
 Function: splitOnProcs

 Parameters:
 symbmtx    - Symbolic matrix
 extracblk  -
 extracost  -
 ctrl       -
 dofptr     -
 cblknum    -
 procnbr    -
 */
void
splitSmart( const BlendCtrl    *ctrl,
            const SymbolMatrix *symbmtx,
            ExtraCblk_t        *extracblk,
            Cand               *candtab)
{
    SymbolBlok   *curblok;
    pastix_int_t *nblocksperline = NULL;
    pastix_int_t  cblknum, bloknum, line;
    pastix_int_t  fsplitrow = -1;
    pastix_int_t  method, authorized_percent;
#define SPLITSYMBOL_METHOD_DEFAULT  0
#define SPLITSYMBOL_METHOD_CONSTANT 1
#define SPLITSYMBOL_METHOD_UPPER    2
    method = SPLITSYMBOL_METHOD_DEFAULT;
    if (pastix_blend_with_constant_split())
        method = SPLITSYMBOL_METHOD_CONSTANT;
    else if (pastix_blend_with_smallest_upper_split())
        method = SPLITSYMBOL_METHOD_UPPER;
    authorized_percent = pastix_blend_split_percent();

    for(cblknum = 0; cblknum<symbmtx->cblknbr; cblknum++)
    {
        pastix_int_t fcolnum = symbmtx->cblktab[cblknum].fcolnum;
        pastix_int_t lcolnum = symbmtx->cblktab[cblknum].lcolnum;
        pastix_int_t candnbr;
        pastix_int_t step;
        pastix_int_t nseq;
        pastix_int_t width;

        /*
         * Compute the number of cblk to be generated by split,
         * for instance we choose to split at the maximum
         */
        candnbr = candtab[ cblknum ].lcandnum
            -     candtab[ cblknum ].fcandnum + 1;

        width = lcolnum - fcolnum + 1;

        nseq = computeNbSplit( ctrl, candnbr, width );
        if (nseq <= 1)
            continue;

        if ( fsplitrow == -1 ) {
            fsplitrow = fcolnum;
            nblocksperline = computeNbBlocksPerLine( symbmtx, fsplitrow );
            nblocksperline -= fsplitrow;
        }

        /* Adapt the step to the segments number */
        step = pastix_iceil( width,  nseq );
        assert( step > 0 );
        nseq--;

        /* { */
        /*     pastix_int_t t, tolerance = 0, min = symbmtx->bloknbr; */

        /*     for(t = symbmtx->cblktab[cblknum].fcolnum; */
        /*         t < symbmtx->cblktab[cblknum].lcolnum; t++) */
        /*     { */
        /*         tolerance += nblocksperline[t]; */
        /*         min = pastix_imin( min, nblocksperline[t] ); */
        /*     } */
        /*     tolerance /= (width-1); */
        /*     pastix_print( 0, 0, "Split %-5ld: Split min=%ld, avg=%ld, width=%ld, nseq=%ld, step=%ld: (", */
        /*                   cblknum, min, tolerance, width, nseq+1, step ); */
        /* } */

        /* Create the new cblk */
        {
            pastix_int_t fcol, lcol;
            pastix_int_t nbcblk = 0;

            fcol = fcolnum;
            while( fcol <= lcolnum )
            {
                if (SPLITSYMBOL_METHOD_CONSTANT == method) {
                    lcol = fcol + computeConstantSplit( nblocksperline + fcol,
                                                        step, width,
                                                        authorized_percent );
                }
                else if (SPLITSYMBOL_METHOD_UPPER == method) {
                    lcol = fcol + computeSmallestSplit_max( nblocksperline + fcol,
                                                            step, width,
                                                            authorized_percent );
                }
                else {
                    lcol = fcol + computeSmallestSplit( nblocksperline + fcol,
                                                        step, width,
                                                        authorized_percent );
                }

                assert( (lcol > fcol) && (lcol <= lcolnum) );

                extraCblkAdd( extracblk, fcol, lcol );
                nbcblk++;

                /* pastix_print( 0, 0, "(%ld, %ld) ", */
                /*               nblocksperline[lcol], (lcol-fcol+1) ); */

                width = width - (lcol - fcol + 1);
                fcol = lcol + 1;
            }

            /*
             * Mark the cblk as being splitted
             */
            extracblk->addcblk += nbcblk-1;
            extracblk->sptcblk[cblknum] = extracblk->curcblk - nbcblk + 1;
            extracblk->sptcbnb[cblknum] = nbcblk;

            /* Update the number of blocks per line*/
            curblok = &(symbmtx->bloktab[symbmtx->cblktab[cblknum].bloknum + 1]) ;
            for(bloknum = symbmtx->cblktab[cblknum].bloknum + 1;
                bloknum < symbmtx->cblktab[cblknum+1].bloknum; bloknum++, curblok++)
            {
                for(line = curblok->frownum; line < curblok->lrownum; line++ )
                {
                    nblocksperline[ line ] += nbcblk-1;
                }
            }
        }
        /* pastix_print( 0, 0, ") \n" ); */
    }

    if ( fsplitrow != -1) {
        nblocksperline += fsplitrow;
        memFree_null( nblocksperline );
    }
}


/*
 Function: splitSymbol

 Repartitioning of the initial symbolic factorization
 and processing of candidate processors group for
 each colum bloc

 Parameters:
 symbmtx - Symbolic matrix.
 ctrl    -
 dofptr  -
 */
void splitSymbol( BlendCtrl    *ctrl,
                  SymbolMatrix *symbmtx )
{
    ExtraCblk_t extracblk;

    /* Init structure to store extra cblks */
    extraCblkInit( symbmtx->cblknbr, &extracblk );

    splitSmart( ctrl, symbmtx, &extracblk, ctrl->candtab );

    /* Merge the initial matrix and the newly generated cblks */
    extraCblkMerge( &extracblk, symbmtx, &(ctrl->candtab) );
    extraCblkExit(&extracblk);

    /* Check that the generated symbol matrix is correct */
    if (ctrl->debug)
        symbolCheck(symbmtx);

    if ( ctrl->clustnum == 0 ) {
        if (ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            symbolPrintStats( symbmtx );
    }

    /* Rk: addcblk field is not erased by Exit call, so we can freely use it */
    if ( extracblk.addcblk )
    {
        /* Update cost matrix to fill-in blank of newly generated blocks */
        costMatrixExit(ctrl->costmtx);
        memFree_null(ctrl->costmtx);
        ctrl->costmtx = costMatrixBuild( symbmtx,
                                         ctrl->iparm[IPARM_FLOAT],
                                         ctrl->iparm[IPARM_FACTORIZATION] );

        if (ctrl->updatecandtab)
        {
            /* Update elimination tree */
            if (ctrl->etree != NULL)
                eTreeExit(ctrl->etree);

            ctrl->etree = eTreeBuild(symbmtx);

            /* Initialize costs in elimination tree and candtab array for proportionnal mapping */
            candBuild( ctrl->autolevel,
                       ctrl->level2D,
                       ctrl->ratiolimit,
                       ctrl->candtab,
                       ctrl->etree,
                       symbmtx,
                       ctrl->costmtx );
        }
    }
}
