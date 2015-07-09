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
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "smart_cblk_split.h"

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
                                    5);
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

        nseq = width / blas_max_col; //pastix_iceil( width, blas_max_col );
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
            step = width / (abs * candnbr); //pastix_iceil( width, (abs * candnbr) );

            step = pastix_imax(step, blas_min_col);
            step = pastix_imin(step, blas_max_col);

            /* Ceil */
            nseq = width / step; //pastix_iceil( width, step );
        }
        else
        {
            nseq = pastix_iceil( width, blas_min_col );
        }
    }

    /* Make sure cblk are at least blas_min_col wide */
    if ( nseq > 1 && (width / nseq) < blas_min_col ) {
        nseq--;
    }

    return nseq;
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
void splitOnProcs2( const BlendCtrl    *ctrl,
                    const SymbolMatrix *symbmtx,
                    ExtraCblk_t        *extracblk,
                    Cand               *candtab)
{
    pastix_int_t i, cblknum;

    for(cblknum = 0; cblknum<symbmtx->cblknbr; cblknum++)
    {
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

        width = symbmtx->cblktab[cblknum].lcolnum
            -   symbmtx->cblktab[cblknum].fcolnum + 1;

        nseq = computeNbSplit( ctrl, candnbr, width );
        if (nseq < 4)
            continue;

        /* Adapt the step to the segments number */
        step = width / nseq; //pastix_iceil( width,  nseq );
        assert( step > 0 );
        nseq--;

        /* Create the new cblk */
#ifdef SMART_CBLK_SPLIT
        {
            pastix_int_t fcolnum = symbmtx->cblktab[cblknum].fcolnum;
            pastix_int_t *seq;
            smart_cblk_split(ctrl,
                             symbmtx,
                             cblknum,
                             candnbr,
                             ctrl->blcolmin,
                             ctrl->blcolmax,
                             &nseq,
                             &seq);

            for(i=0;i<nseq;i++)
            {
                extraCblkAdd( extracblk,
                              fcolnum + seq[2*i],
                              fcolnum + seq[2*i+1] );
            }

            extraCblkAdd( extracblk,
                          fcolnum + seq[2*nseq],
                          symbmtx->cblktab[cblknum].lcolnum );
        }
#else
        {
            pastix_int_t fcolnum = symbmtx->cblktab[cblknum].fcolnum;

            for(i=0;i<nseq;i++)
            {
                extraCblkAdd( extracblk,
                              fcolnum + step * i,
                              fcolnum + step * (i+1) -1 );
            }

            extraCblkAdd( extracblk,
                          fcolnum + step * nseq,
                          symbmtx->cblktab[cblknum].lcolnum );
        }
#endif
        /*
         * Mark the cblk as being splitted
         */
        extracblk->addcblk += nseq;
        extracblk->sptcblk[cblknum] = extracblk->curcblk - nseq;
        extracblk->sptcbnb[cblknum] = nseq+1;
    }
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
        /*     pastix_print( 0, 0, "Tolerance of cblk %ld is %ld and min is %ld: (", */
        /*                   cblknum, tolerance, min ); */
        /* } */

        /* Create the new cblk */
        {
            pastix_int_t fcol, lcol;
            pastix_int_t nbcblk = 0;

            fcol = fcolnum;
            while( fcol <= lcolnum )
            {
                if (SPLITSYMBOL_METHOD_CONSTANT == method)
                    lcol = fcol + computeConstantSplit( nblocksperline + fcol,
                                                        step, width,
                                                        authorized_percent );
                else if (SPLITSYMBOL_METHOD_UPPER == method)
                    lcol = fcol + computeSmallestSplit_max( nblocksperline + fcol,
                                                            step, width,
                                                            authorized_percent );
                else
                    lcol = fcol + computeSmallestSplit( nblocksperline + fcol,
                                                        step, width,
                                                        authorized_percent );

                assert( (lcol > fcol) && (lcol <= lcolnum) );

                extraCblkAdd( extracblk, fcol, lcol );
                nbcblk++;

                /* pastix_print( 0, 0, "(%ld, %ld, %ld) ", */
                /*               nblocksperline[fcol], (lcol-fcol+1), step ); */

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

    /* Stupid split */
    if (1)
    {
#ifdef SMART_CBLK_SPLIT
        /* Build the elimination graph from the new symbolic partition */
        {
            Clock timer_current;
            pastix_print(0 , 0, OUT_BLEND_ELIMGRAPH2 );
            clockStart(timer_current);

            MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);
            eGraphInit(ctrl->egraph);
            eGraphBuild(ctrl->egraph, symbmtx);

            clockStop(timer_current);
            pastix_print( 0, 0, "--Graph build at time: %g --\n", clockVal(timer_current) );
        }
#endif
        splitOnProcs2( ctrl, symbmtx, &extracblk, ctrl->candtab );
#ifdef SMART_CBLK_SPLIT
        eGraphExit(ctrl->egraph);
#endif
    }
    else
    {
        splitSmart( ctrl, symbmtx, &extracblk, ctrl->candtab );
    }

    /* Merge the initial matrix and the newly generated cblks */
    extraCblkMerge( &extracblk, symbmtx, &(ctrl->candtab) );

    /* Check that the generated symbol matrix is correct */
    if (ctrl->debug)
        symbolCheck(symbmtx);


    if (ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    {
        pastix_int_t i, j;
        double block_height_sum = 0.0;
        double cblk_width_sum = 0.0;
        pastix_int_t block_height_min = (symbmtx->cblktab[0].lcolnum - symbmtx->cblktab[0].fcolnum + 1);
        pastix_int_t cblk_width_min = (symbmtx->bloktab[0].lrownum - symbmtx->bloktab[0].frownum + 1);;
        pastix_int_t block_height_max = 0.0;
        pastix_int_t cblk_width_max = 0.0;

        pastix_print( ctrl->clustnum, 0,
                      "Number of column blocks modified by splitting   : %ld\n"
                      "Number of column blocks created by splitting    : %ld\n"
                      "Number of blocks creating by splitting          : %ld\n"
                      "Number of blocks creating by splitting of fcblk : %ld\n"
                      "Oldsymbol cblknbr = %ld, bloknbr = %ld\n"
                      "Newsymbol cblknbr = %ld, bloknbr = %ld\n",
                      (long int)(extracblk.curcblk + 1),
                      (long int)(extracblk.addcblk),
                      (long int)(extracblk.addblok),
                      (long int)(extracblk.addblof),
                      (long int)(symbmtx->cblknbr - extracblk.addcblk),
                      (long int)(symbmtx->bloknbr - extracblk.addblok),
                      (long int)(symbmtx->cblknbr),
                      (long int)(symbmtx->bloknbr) );

        for (j = 0; j < symbmtx->cblknbr; j++)
        {
            pastix_int_t cblk_width = ( symbmtx->cblktab[j].lcolnum -
                                        symbmtx->cblktab[j].fcolnum + 1);
            cblk_width_sum += (double)cblk_width;
            cblk_width_max = pastix_imax(cblk_width_max, cblk_width);
            cblk_width_min = pastix_imin(cblk_width_min, cblk_width);
            for (i = symbmtx->cblktab[j].bloknum+1; i < symbmtx->cblktab[j+1].bloknum; i++)
            {
                pastix_int_t block_height = ( symbmtx->bloktab[i].lrownum -
                                              symbmtx->bloktab[i].frownum + 1);
                block_height_sum += (double)block_height;
                block_height_max = pastix_imax(block_height_max, block_height);
                block_height_min = pastix_imin(block_height_min, block_height);
            }
        }
        fprintf(stdout, "Cblk width : Avg %6.3g\t Min %ld\t Max %ld\n",
                cblk_width_sum/symbmtx->cblknbr,
                cblk_width_min, cblk_width_max);
        fprintf(stdout, "Extra diagonal block height : Avg %6.3g\t Min %ld\t Max %ld\n",
                block_height_sum/(symbmtx->bloknbr-symbmtx->cblknbr),
                block_height_min, block_height_max);
        if (ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
        {
            pastix_int_t cblk_step = (cblk_width_max - cblk_width_min)/10;
            pastix_int_t block_step = (block_height_max - block_height_min)/10;
            pastix_int_t cblk_counters[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            pastix_int_t block_counters[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            if (cblk_step == 0) cblk_step++;
            if (block_step == 0) block_step++;
            if (cblk_width_min + cblk_step * 10 < cblk_width_max)
                cblk_step++;
            if (block_height_min + block_step * 10 < block_height_max)
                block_step++;
            for (j = 0; j < symbmtx->cblknbr; j++)
            {
                pastix_int_t cblk_width = ( symbmtx->cblktab[j].lcolnum -
                                            symbmtx->cblktab[j].fcolnum + 1);
                cblk_counters[(cblk_width-cblk_width_min)/cblk_step]++;
                for (i = symbmtx->cblktab[j].bloknum+1; i < symbmtx->cblktab[j+1].bloknum; i++)
                {
                    pastix_int_t block_height = ( symbmtx->bloktab[i].lrownum -
                                                  symbmtx->bloktab[i].frownum + 1);
                    block_counters[(block_height-block_height_min)/block_step]++;
                }
            }
            fprintf(stdout, "Cblk width range :");
            for (j = 0; j < 10; j ++)
                fprintf(stdout, " [ %3ld - %3ld ]",
                        cblk_width_min+j*cblk_step,
                        cblk_width_min+(j+1)*cblk_step);
            fprintf(stdout, "\n");
            fprintf(stdout, "Number of Cblks :  ");
            for (j = 0; j < 10; j ++)
                fprintf(stdout, "[ %9ld ]",
                        cblk_counters[j]);
            fprintf(stdout, "\n");

            fprintf(stdout, "Block height range :");
            for (j = 0; j < 10; j ++)
                fprintf(stdout, "[ %3ld - %3ld ]",
                        block_height_min+j*block_step,
                        block_height_min+(j+1)*block_step);
            fprintf(stdout, "\n");
            fprintf(stdout, "Number of Blocks :   ");
            for (j = 0; j < 10; j ++)
                fprintf(stdout, "[ %9ld ]",
                        block_counters[j]);

            fprintf(stdout, "\n");
        }
    }
    extraCblkExit(&extracblk);

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
