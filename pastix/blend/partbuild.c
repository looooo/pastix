#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "extrastruct.h"
#include "elimin.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "blendctrl.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "partbuild.h"

#define OLD_MERGE
/*+
 Make a new SymbolMatrix and CostMatrix from the former ones
 and the Extra ones (that contains splitted bloks)
 +*/
/* OIMBE : retravailler partBuild avec des memcpy */

void symbolMerge( BlendCtrl *ctrl, const Dof * dofptr,
                  SymbolMatrix *newsymb, ExtraSymbolMatrix *extrasymb,
                  CostMatrix   *newcost, ExtraCostMatrix   *extracost,
                  Cand        **candtab )
{
    pastix_int_t  i, j, k, s;
    pastix_int_t  curbloknum;
    pastix_int_t  sptbloknum;
    pastix_int_t *newnum      = NULL;
    pastix_int_t *extranewnum = NULL;
    pastix_int_t  facing_splitted_cnt = 0;

    SymbolMatrix *oldsymb;
    CostMatrix   *oldcost;
    Cand         *oldcand = *candtab;
    Cand         *newcand;

    /* No splitted cblk: partition remains the same */
    fprintf( stderr, "Number of generated block is %ld\n", extrasymb->curcblk );
    if( extrasymb->curcblk == 0 ) {
        return;
    }

    pastix_print( ctrl->clustnum, 0,
                  "Number of column blocks created by splitting : %d\n"
                  "Number of blocks creating by splitting       : %d\n",
                  (int)(extrasymb->addcblk), (int)(extrasymb->addblok));

    /* Allocate new symbol */
    MALLOC_INTERN(oldsymb, 1, SymbolMatrix);
    memcpy( oldsymb, newsymb, sizeof(SymbolMatrix) );

    newsymb->cblknbr = oldsymb->cblknbr + extrasymb->addcblk;
    newsymb->bloknbr = oldsymb->bloknbr + extrasymb->addblok;
    MALLOC_INTERN(newsymb->cblktab, newsymb->cblknbr+1, SymbolCblk);
    MALLOC_INTERN(newsymb->bloktab, newsymb->bloknbr,   SymbolBlok);
    memset( newsymb->bloktab, 0, newsymb->bloknbr * sizeof(SymbolBlok) );

    /* Allocate new CostMatrix */
    MALLOC_INTERN(oldcost, 1, CostMatrix);
    memcpy( oldcost, newcost, sizeof(CostMatrix) );

    costInit(newcost);
    MALLOC_INTERN(newcost->bloktab, newsymb->bloknbr, CostBlok);

    /* Allocate new candtab */
    MALLOC_INTERN(newcand, newsymb->cblknbr, Cand);
    memset( newcand, 0, newsymb->cblknbr * sizeof(Cand) );

#if defined(OLD_MERGE)
    /*
     * We use the sptcbnb array to get the new numbering of the former cblk
     * in the new symbolic matrix
     * newnum[i+1] becomes the new number of the first cblk generated from the
     * split of former cblk number i.
     */
    MALLOC_INTERN(newnum, oldsymb->cblknbr, pastix_int_t);
    memcpy(newnum+1, extrasymb->sptcbnb, (oldsymb->cblknbr-1) * sizeof(pastix_int_t));
    newnum[0] = 0;
    for(i=1; i<oldsymb->cblknbr; i++) {
        newnum[i] += newnum[i-1];
        assert( (newnum[i] >= 0) && (newnum[i] < newsymb->cblknbr) );
        assert( newnum[i] > newnum[i-1] );
    }

    /*
     * Now, we use sptcblk and newnum to get the new numbering of all generated
     * cblk owned by the extra symbolic matrix
     */
    MALLOC_INTERN(extranewnum, extrasymb->curcblk, pastix_int_t);
    for(i=0; i<oldsymb->cblknbr; i++)
    {
        if(extrasymb->sptcblk[i] >= 0)
        {
            for(j=0; j<extrasymb->sptcbnb[i]; j++) {
                extranewnum[extrasymb->sptcblk[i]+j] = newnum[i]+j;
                assert( (extranewnum[extrasymb->sptcblk[i]+j] >= 0) &&
                        (extranewnum[extrasymb->sptcblk[i]+j] < newsymb->cblknbr) );
                assert( extranewnum[extrasymb->sptcblk[i]+j] >= newnum[i] );
            }
        }
    }

#else
    /*
     * We use the sptcbnb array to get the new numbering of the former cblk
     * in the new symbolic matrix
     * newnum[i+1] becomes the new number of the first cblk generated from the
     * split of former cblk number i.
     */
    MALLOC_INTERN(newnum,      oldsymb->cblknbr,   pastix_int_t);
    MALLOC_INTERN(extranewnum, extrasymb->curcblk, pastix_int_t);
    memcpy(newnum+1, extrasymb->sptcbnb, (oldsymb->cblknbr-1) * sizeof(pastix_int_t));
    {
        pastix_int_t sptcbnb, sptcblk;
        pastix_int_t nbcblk2copy = 0;
        pastix_int_t firstcblk   = 0;

        /* First cblk */
        newnum[0] = 0;
        sptcbnb = extrasymb->sptcbnb[0];
        assert( sptcbnb > 0 );

        /* Splitted cblk */
        if (sptcbnb > 1)
        {
            sptcblk = extrasymb->sptcblk[0];
            assert( sptcblk >= 0 );

            for(j=0; j<sptcbnb; j++, sptcblk++) {
                extranewnum[sptcblk] = j;

                assert( (extranewnum[sptcblk] >= 0) &&
                        (extranewnum[sptcblk] < newsymb->cblknbr) );

                memcpy( newcand + extranewnum[ sptcblk ],
                        oldcand, sizeof(Cand) );
            }
        }
        else {
            nbcblk2copy++;
        }

        for(i=1; i<oldsymb->cblknbr; i++) {

            sptcbnb = extrasymb->sptcbnb[i];
            newnum[i] += newnum[i-1];

            assert( sptcbnb > 0 );
            assert( (newnum[i] >= 0) && (newnum[i] < newsymb->cblknbr) );

            /* Splitted cblk */
            if (sptcbnb > 1)
            {
                /* Copy the unchanged cblk */
                if ( nbcblk2copy > 0 ) {

                    memcpy( newsymb->cblktab + newnum[ firstcblk ],
                            oldsymb->cblktab + firstcblk,
                            nbcblk2copy * sizeof(SymbolCblk) );

                    memcpy( newcand + newnum[ firstcblk ],
                            oldcand + firstcblk,
                            nbcblk2copy * sizeof(Cand) );

                    nbcblk2copy = 0;
                }

                sptcblk = extrasymb->sptcblk[i];
                assert( sptcblk >= 0 );

                for(j=0; j<sptcbnb; j++, sptcblk++) {
                    extranewnum[sptcblk] = newnum[i]+j;

                    assert( (extranewnum[sptcblk] >= 0) &&
                            (extranewnum[sptcblk] <  newsymb->cblknbr) );
                    assert(  extranewnum[sptcblk] >= newnum[i] );

                    memcpy( newcand + extranewnum[ sptcblk ],
                            oldcand + i, sizeof(Cand) );
                }
            }
            else
            {
                nbcblk2copy++;
            }
        }
    }
    memFree_null(oldcand);
#endif

    /* Fill in the new symbolic matrix resulting from the splitting of the former one */
    curbloknum = 0;
    for(i=0; i<oldsymb->cblknbr; i++)
    {
        pastix_int_t newcblknum;

        if(extrasymb->sptcblk[i] < 0) /* not a splitted cblk */
        {
            newcblknum = newnum[i];
            assert( newcblknum < newsymb->cblknbr );

#if defined(OLD_MERGE)
	    memcpy( newcand + newcblknum, oldcand + i, sizeof(Cand) );
	    newsymb->cblktab[newcblknum].fcolnum = oldsymb->cblktab[i].fcolnum;
	    newsymb->cblktab[newcblknum].lcolnum = oldsymb->cblktab[i].lcolnum;
#endif
            /* Update bloknum */
            newsymb->cblktab[newcblknum].bloknum = curbloknum;

            for( j = oldsymb->cblktab[i].bloknum;
                 j < oldsymb->cblktab[i+1].bloknum; j++ )
            {
                /*
                 * Not a splitted blok, so its facing diag is not splitted
                 * NB: even if a blok is not splitted while its facing cblk is
                 * splitted , it's considered as splitted
                 */
                if(extrasymb->sptblok[j] < 0)
                {
                    pastix_int_t fcblknum = oldsymb->bloktab[j].cblknum;
                    assert( fcblknum   < oldsymb->cblknbr );
                    fcblknum = newnum[ fcblknum ];

                    newsymb->bloktab[curbloknum].frownum = oldsymb->bloktab[j].frownum;
                    newsymb->bloktab[curbloknum].lrownum = oldsymb->bloktab[j].lrownum;
                    newsymb->bloktab[curbloknum].cblknum = fcblknum;
                    newsymb->bloktab[curbloknum].levfval = oldsymb->bloktab[j].levfval;

                    assert( (fcblknum >= 0) &&
                            (fcblknum < newsymb->cblknbr) );

                    newcost->bloktab[curbloknum].contrib = oldcost->bloktab[j].contrib;
                    newcost->bloktab[curbloknum].linenbr = oldcost->bloktab[j].linenbr;
                    curbloknum++;
		    assert( curbloknum <= newsymb->bloknbr );
                }
                /* Splitted blok in a non splitted cblk -> the facing diagblok is splitted */
                else
                {
                    facing_splitted_cnt += extrasymb->sptblnb[j]-1;
                    for(k = extrasymb->sptblok[j];
                        k < extrasymb->sptblok[j] + extrasymb->sptblnb[j]; k++)
                    {
                        pastix_int_t fcblknum = extranewnum[ extrasymb->bloktab[k].cblknum ];

			assert( curbloknum < newsymb->bloknbr );
                        newsymb->bloktab[curbloknum].frownum = extrasymb->bloktab[k].frownum;
                        newsymb->bloktab[curbloknum].lrownum = extrasymb->bloktab[k].lrownum;
                        newsymb->bloktab[curbloknum].cblknum = fcblknum;
                        newsymb->bloktab[curbloknum].levfval = extrasymb->bloktab[k].levfval;

                        assert( (fcblknum >= 0) &&
                                (fcblknum < newsymb->cblknbr) );

                        newcost->bloktab[curbloknum].contrib = extracost->bloktab[k].contrib;
                        newcost->bloktab[curbloknum].linenbr = extracost->bloktab[k].linenbr;

                        curbloknum++;
			assert( curbloknum <= newsymb->bloknbr );
                    }
                }
            }
        }
        /* Splitted cblk */
        else
        {
            for(j = extrasymb->sptcblk[i];
                j < extrasymb->sptcblk[i] + extrasymb->sptcbnb[i]; j++)
            {
                newcblknum = extranewnum[j];
                assert( newcblknum < newsymb->cblknbr );

#if defined(OLD_MERGE)
		memcpy( newcand + newcblknum, oldcand + i, sizeof(Cand) );
#endif

                newsymb->cblktab[newcblknum].fcolnum = extrasymb->cblktab[j].fcolnum;
                newsymb->cblktab[newcblknum].lcolnum = extrasymb->cblktab[j].lcolnum;
                newsymb->cblktab[newcblknum].bloknum = curbloknum;

                /* Treat blok created by splitting of the diag blok */
                for(k = extrasymb->cblktab[j].bloknum;
                    k <(extrasymb->cblktab[j].bloknum + extrasymb->sptcblk[i] + extrasymb->sptcbnb[i] - j); k++)
                {
                    pastix_int_t fcblknum = extrasymb->bloktab[k].cblknum;

		    assert( curbloknum < newsymb->bloknbr );
                    newsymb->bloktab[curbloknum].frownum = extrasymb->bloktab[k].frownum;
                    newsymb->bloktab[curbloknum].lrownum = extrasymb->bloktab[k].lrownum;
                    newsymb->bloktab[curbloknum].cblknum = extranewnum[fcblknum];
                    newsymb->bloktab[curbloknum].levfval = extrasymb->bloktab[k].levfval;

                    assert(newsymb->bloktab[curbloknum].frownum >= extrasymb->cblktab[fcblknum].fcolnum);
                    assert(newsymb->bloktab[curbloknum].lrownum <= extrasymb->cblktab[fcblknum].lcolnum);

                    curbloknum++;
		    assert( curbloknum <= newsymb->bloknbr );
                }

                sptbloknum = k;
                for(k = oldsymb->cblktab[i].bloknum+1;
                    k < oldsymb->cblktab[i+1].bloknum; k++)
                {
                    if(extrasymb->sptblok[k] < 0)
                    {
                        pastix_int_t fcblknum = extrasymb->bloktab[sptbloknum].cblknum;
                        assert( fcblknum < oldsymb->cblknbr );
                        fcblknum = newnum[ fcblknum ];

			assert( curbloknum < newsymb->bloknbr );
                        newsymb->bloktab[curbloknum].frownum = oldsymb->bloktab[k].frownum;
                        newsymb->bloktab[curbloknum].lrownum = oldsymb->bloktab[k].lrownum;
                        newsymb->bloktab[curbloknum].cblknum = fcblknum;
                        newsymb->bloktab[curbloknum].levfval = oldsymb->bloktab[k].levfval;
                        sptbloknum++;
                        curbloknum++;
			assert( curbloknum <= newsymb->bloknbr );
                    }
                    else
                    {
                        facing_splitted_cnt += extrasymb->sptblnb[k]-1;
                        for(s=extrasymb->sptblok[k];s<extrasymb->sptblok[k]+extrasymb->sptblnb[k];s++)
                        {
                            pastix_int_t fcblknum = extranewnum[extrasymb->bloktab[s].cblknum];

			    assert( curbloknum < newsymb->bloknbr );
                            newsymb->bloktab[curbloknum].frownum = extrasymb->bloktab[s].frownum;
                            newsymb->bloktab[curbloknum].lrownum = extrasymb->bloktab[s].lrownum;
                            newsymb->bloktab[curbloknum].cblknum = fcblknum;
                            newsymb->bloktab[curbloknum].levfval = extrasymb->bloktab[s].levfval;
                            sptbloknum++;
                            curbloknum++;
			    assert( curbloknum <= newsymb->bloknbr );
                        }
                    }
                }
            }
        }
    }

    assert(curbloknum == newsymb->bloknbr);

    /* Free old versions and temporary buffer */
    symbolExit(oldsymb);
    memFree_null(oldsymb);
    costExit(oldcost);

#if defined(OLD_MERGE)
    memFree_null(oldcand);
#endif
    memFree_null(newnum);
    memFree_null(extranewnum);

    /* Virtual cblk to avoid side effect in the loops on cblk bloks */
    newsymb->cblktab[newsymb->cblknbr].fcolnum = newsymb->cblktab[newsymb->cblknbr-1].lcolnum+1;
    newsymb->cblktab[newsymb->cblknbr].lcolnum = newsymb->cblktab[newsymb->cblknbr-1].lcolnum+1;
    newsymb->cblktab[newsymb->cblknbr].bloknum = curbloknum;

    /* TODO: update on the fly */
    /* We have to compute the cost of a splitted cblk */
    for(i=0; i<newsymb->cblknbr; i++)
    {
        cblkComputeCost(i, newcost, newsymb, dofptr);
    }

    if (ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    {
        double        block_height_sum = 0.0;
        double        cblk_width_sum = 0.0;
        for (j = 0; j < newsymb->cblknbr; j++)
        {
            cblk_width_sum += (double)(newsymb->cblktab[j].lcolnum - newsymb->cblktab[j].fcolnum + 1);

            for (i = newsymb->cblktab[j].bloknum+1; i < newsymb->cblktab[j+1].bloknum; i++)
            {
                block_height_sum += (double)(newsymb->bloktab[i].lrownum - newsymb->bloktab[i].frownum + 1);
            }
        }
        fprintf(stdout, "Average cblk size : %g\n", cblk_width_sum/newsymb->cblknbr);
        fprintf(stdout, "Average extra diagonal block height : %g\n", block_height_sum/(newsymb->bloknbr-newsymb->cblknbr));
        fprintf(stdout, "Number of blocks created due to facing block splitting : %d\n", (int)facing_splitted_cnt);
    }

    symbolCheck(newsymb);

    *candtab = newcand;
    return;
}
