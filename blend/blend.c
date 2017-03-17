/*
  File: blend.c

  Main blend source code.

  Re-split column blocs and distribute them on processors.

  Authors:
    - Pascal  Henon
    - Pierre  Ramet
    - Mathieu Faverge
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "out.h"
#include "ftgt.h"
#include "cost.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "elimin.h"
#include "extendVector.h"
#include "cand.h"
#include "blendctrl.h"
#include "simu.h"
/* #include "dof.h" */
/* #include "costfunc.h" */
/* #include "splitpartlocal.h" */
#include "solver_check.h"
#include "task.h"
#include "solver_check.h"
#include "fanboth2.h"
#include "blend.h"
#include "order.h"

/* #include "assert.h" */

/*
 * Function: solverBlend
 *
 * Main blend function
 *
 * Build the elimination graph from the symbolic partition.
 *
 * Build the cost matrix from the symbolic partition.
 *
 * Build the elimination tree from the symbolic partition.
 *
 * Distribute each column bloc on candidate processors.
 *
 * Build a new symbol matrix...
 *
 * Parameters:
 *   solvmtx    - Solver matrix structure.
 */
void solverBlend(BlendCtrl    *ctrl,
                 SolverMatrix *solvmtx,
                 SymbolMatrix *symbmtx)
{
    SimuCtrl     *simuctrl;
    double        timer_all     = 0.;
    double        timer_current = 0.;
    pastix_int_t  clustnum = ctrl->clustnum;
    pastix_int_t  clustnbr = ctrl->clustnbr;

    clockStart(timer_all);

    if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print( clustnum, 0, OUT_BLEND_CONF,
                      (long)clustnbr, (long)ctrl->local_nbcores, (long)ctrl->local_nbthrds);

    /* Verify the coherence of the initial symbol matrix */
    if(ctrl->debug)
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_CHKSMBMTX );
        }
        symbolCheck(symbmtx);
    }

    /* Build the elimination tree from the symbolic partition */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES) {
            pastix_print( clustnum, 0, OUT_BLEND_ELIMTREE );
        }
        clockStart(timer_current);

        ctrl->etree = eTreeBuild(symbmtx);

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_ELIMTREE_TIME,
                          clockVal(timer_current) );
        }
    }

    /* Build the cost matrix from the symbolic partition */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_COSTMATRIX );
        }
        clockStart(timer_current);

        ctrl->costmtx = costMatrixBuild( symbmtx,
                                         ctrl->iparm[IPARM_FLOAT],
                                         ctrl->iparm[IPARM_FACTORIZATION] );

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_COSTMATRIX_TIME,
                          clockVal(timer_current) );
        }
    }

    /* Build the candtab array to store candidate information on each cblk */
    {
        MALLOC_INTERN(ctrl->candtab, symbmtx->cblknbr, Cand);
        candInit( ctrl->candtab, symbmtx->cblknbr );

        /* Initialize costs in elimination tree and candtab array for proportionnal mapping */
        candBuild( ctrl->autolevel,
                   ctrl->level2D,
                   ctrl->ratiolimit,
                   ctrl->candtab,
                   ctrl->etree,
                   symbmtx,
                   ctrl->costmtx );
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_ELIMTREE_TOTAL_COST,
                          ctrl->etree->nodetab[ eTreeRoot(ctrl->etree) ].subtree );
        }
    }

    /* Proportional mapping step that distributes the candidates over the tree */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_PROPMAP );
        }
        clockStart(timer_current);

        propMappTree( ctrl->candtab,
                      ctrl->etree,
                      ctrl->total_nbcores,
                      ctrl->nocrossproc,
                      ctrl->allcand );

        /* Set the cluster candidates according to the processor candidates */
        candSetClusterCand( ctrl->candtab, symbmtx->cblknbr,
                            ctrl->core2clust, ctrl->total_nbcores );

        clockStop(timer_current);

        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_PROPMAP_TIME,
                          clockVal(timer_current) );
        }

        /* Let's check the result if ask */
        if ( ctrl->debug ) {
            assert( candCheck( ctrl->candtab, symbmtx ) );
        }
    }

    /*
     * Split the existing symbol matrix according to the number of candidates
     * and cblk types.
     * It takes the original symbol and candtab, and return the new symbol and
     * candtab. If the symbmtx is modified, the costmtx is updated, as well as
     * the tree.
     */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_SPLITSYMB );
        }
        clockStart(timer_current);

        splitSymbol(ctrl, symbmtx);

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_SPLITSYMB_TIME,
                          clockVal(timer_current) );
        }
    }

    if(ctrl->count_ops && (ctrl->leader == clustnum)) {
        double thflops = 0.;
        symbolGetFlops( symbmtx,
                        ctrl->iparm[IPARM_FLOAT],
                        ctrl->iparm[IPARM_FACTORIZATION],
                        &thflops, &(ctrl->dparm[DPARM_FACT_RLFLOPS]) );
        // TODO: check why the splitsymbol changes the number of flops, it should not.
        //assert( thflops == ctrl->dparm[DPARM_FACT_THFLOPS] );
    }

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbolblend.eps", "w");
        symbolDraw(symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    if (0)
    {
        FILE *file = fopen("symbgen2", "w");
        symbolSave( symbmtx, file );
        fclose(file);
    }

    /* Simulation step to perform the data distribution over the nodes and compute the priorities of each task */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_BUILDSIMU );
        }
        clockStart(timer_current);

        /* Initialize simulation structure */
        MALLOC_INTERN(simuctrl, 1, SimuCtrl);
        simuInit( simuctrl, symbmtx, ctrl->candtab,
                  ctrl->clustnbr,
                  ctrl->total_nbcores );

        /* Create task array */
        taskBuild(simuctrl, symbmtx, ctrl->candtab);
        clockStop(timer_current);

        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_BUILDSIMU_TIME,
                          clockVal(timer_current),
                          (long)simuctrl->tasknbr );
        }

        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_SIMU );
        }
        clockStart(timer_current);

        simuRun( simuctrl, ctrl, symbmtx );

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_SIMU_TIME,
                          clockVal(timer_current) );
        }
    }

    /* Build the elimination graph from the new symbolic partition */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_ELIMGRAPH );
        }
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);
        eGraphInit(ctrl->egraph);
        eGraphBuild(ctrl->egraph, symbmtx);

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_ELIMGRAPH_TIME,
                          clockVal(timer_current) );
        }
    }

#ifdef PASTIX_DYNSCHED
    /**
     * If dynamic scheduling is asked, let's perform a second proportionnal
     * mapping step:
     *    - this is made only on local data
     *    - no crossing is allowed between branches
     */
    {
        clockStart(timer_current);

        splitPartLocal( ctrl, simuctrl, symbmtx );

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
            pastix_print( clustnum, 0, "  -- Split build at time: %g --\n", clockVal(timer_current));
    }
#endif

    /* CostMatrix and Elimination Tree are no further used */
    costMatrixExit( ctrl->costmtx );
    memFree_null( ctrl->costmtx );
    eTreeExit( ctrl->etree );

    /**
     * Generate the final solver structure that collects data from the different
     * simulation structures and convert to local numbering
     */
    {
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            pastix_print( clustnum, 0, OUT_BLEND_SOLVER );
        }
        clockStart(timer_current);

        solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl);

        clockStop(timer_current);
        if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( clustnum, 0, OUT_BLEND_SOLVER_TIME,
                          clockVal(timer_current) );
            if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
                solverPrintStats( solvmtx );
            }
        }
    }

    /* Free allocated memory */
    simuExit(simuctrl, ctrl->clustnbr, ctrl->total_nbcores, ctrl->local_nbctxts);
    eGraphExit(ctrl->egraph);

    /* Realloc solver memory in a contiguous way  */
    {
        solverRealloc(solvmtx);

#if defined(PASTIX_DEBUG_BLEND)
        if (!ctrl->ricar) {
            if( ctrl->iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
                pastix_print( clustnum, 0, OUT_BLEND_CHKSOLVER );
            }
            solverCheck(solvmtx);
        }
#endif
    }

    /* End timing */
    clockStop(timer_all);
    set_dparm(ctrl->dparm, DPARM_ANALYZE_TIME, clockVal(timer_all) );
}
