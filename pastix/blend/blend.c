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
#include "dof.h"
#include "ftgt.h"
#include "cost.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "solverRealloc.h"
#include "elimin.h"
#include "extrastruct.h"
#include "extendVector.h"
#include "cand.h"
#include "blendctrl.h"
#include "splitpart.h"
#include "simu.h"
#include "costfunc.h"
#include "splitpartlocal.h"
#include "distribPart.h"
#include "solverMatrixGen.h"
#include "solver_check.h"
#include "symbol_cost.h"
#include "task.h"
#include "solver_check.h"
#include "fanboth2.h"
#include "blend.h"
#include "order.h"
#include "bordi.h"

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
                 SymbolMatrix *symbmtx,
                 const Dof    *dofptr)
{
    SimuCtrl     *simuctrl;
    double        timer_all     = 0.;
    double        timer_current = 0.;
    pastix_int_t *bcofind       = NULL;
    pastix_int_t  clustnum = ctrl->clustnum;
    pastix_int_t  clustnbr = ctrl->clustnbr;

    clockStart(timer_all);

    pastix_print( clustnum, 0,
                  OUT_CLUSTNBR "" OUT_PROCNBR "" OUT_THRDNBR,
                  (long)clustnbr, (long)ctrl->local_nbcores, (long)ctrl->local_nbthrds);

    /* Verify the coherence of the initial symbol matrix */
    if(ctrl->debug)
    {
        pastix_print( clustnum, 0, OUT_BLEND_CHKSMBMTX );
        symbolCheck(symbmtx);
    }

    if(ctrl->count_ops && ctrl->leader == clustnum)
        symbCost(ctrl->iparm, ctrl->dparm, symbmtx, dofptr);

    /* Build the elimination tree from the symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_ELIMTREE );
        clockStart(timer_current);

        ctrl->etree = eTreeBuild(symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Tree build at time: %g --\n", clockVal(timer_current));
    }

    /* Build the cost matrix from the symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_COSTMATRIX );
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->costmtx, 1, CostMatrix);
        costMatrixInit(ctrl->costmtx);
        costMatrixBuild(ctrl->costmtx, symbmtx, dofptr);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "-- Cost Matrix build at time:          %g\n",
                      clockVal(timer_current));
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

        pastix_print( clustnum, 0, "-- Total cost of the elimination tree: %g\n",
                      ctrl->etree->nodetab[ eTreeRoot(ctrl->etree) ].subtree);
    }

    /* Proportional mapping step that distributes the candidates over the tree */
    {
        clockStart(timer_current);
        pastix_print( clustnum, 0, "-- Start proportionnal mapping step\n" );

        propMappTree( ctrl->candtab,
                      ctrl->etree,
                      symbmtx,
                      dofptr,
                      ctrl->total_nbcores,
                      ctrl->nocrossproc,
                      ctrl->allcand );

        /* Set the cluster candidates according to the processor candidats */
        candSetClusterCand( ctrl->candtab, symbmtx->cblknbr,
                            ctrl->core2clust, ctrl->total_nbcores );

        /* Let's check the result if ask */
        if (ctrl->debug) {
            assert( candCheck( ctrl->candtab, symbmtx ) );
        }

        clockStop(timer_current);
        pastix_print( clustnum, 0, "-- Proportionnal mapping time:      %g\n", clockVal(timer_current));
    }

    /*
     * Split the existing symbol matrix according to the number of candidates
     * and cblk types.
     * It takes the original symbol and candtab, and return the new symbol and
     * candtab. If the symbmtx is modified, the costmtx is updated, as well as
     * the tree.
     */
    {
        clockStart(timer_current);
        pastix_print( clustnum, 0, "-- Spliting initial partition Version 2\n" );

        splitPart2(ctrl, symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "-- Split build at time: %g --\n", clockVal(timer_current));
    }

    if(ctrl->count_ops && (ctrl->leader == clustnum))
        symbCost(ctrl->iparm, ctrl->dparm, symbmtx, dofptr);

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbolblend.eps", "w");
        symbolDraw(symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    /* Build the elimination graph from the new symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_ELIMGRAPH2 );
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);
        eGraphInit(ctrl->egraph);
        eGraphBuild(ctrl->egraph, symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Graph build at time: %g --\n", clockVal(timer_current) );
    }

    if(ctrl->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    /* initialize simu structure control */
    MALLOC_INTERN(simuctrl, 1, SimuCtrl);
    simuInit(simuctrl, symbmtx,
             ctrl->clustnbr,
             ctrl->total_nbcores,
             symbmtx->cblknbr, symbmtx->bloknbr, ctrl->candtab);

    /* Build tasks */
    if( ctrl->leader == clustnum &&
        ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_TASKGRAPH);
    taskBuild(simuctrl, symbmtx, ctrl->candtab, dofptr, ctrl);

    if( ctrl->leader == clustnum &&
        ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_NBTASK, (long)simuctrl->tasknbr);
    if(ctrl->timer)
    {
        clockStop(timer_current);
        printf("--Task built at time: %g --\n", clockVal(timer_current));
    }

    /* Distribution Phase */
    if(ctrl->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

#ifdef DEBUG_BLEND
    ASSERT(check_candidat(symbmtx, ctrl)>=0,MOD_BLEND);
#endif

    if((ctrl->leader == clustnum) &&
       (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
        fprintf(stdout, OUT_BLEND_DISTPART);
    distribPart(symbmtx, simuctrl, ctrl, dofptr);

    if(ctrl->timer)
    {
        clockStop(timer_current);
        printf("--Distribution computed at time: %g --\n", clockVal(timer_current));
    }

#ifdef PASTIX_DYNSCHED /* 2 eme passe de splitpart */

    if(ctrl->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    /** repartitioning of the initial symbolic factorization
     and processing of candidate processors group for
     each colum bloc **/

    splitPartLocal(ctrl, simuctrl, symbmtx, dofptr);

    if(ctrl->timer)
    {
        clockStop(timer_current);
        printf("--Split build at time: %g --\n", clockVal(timer_current));
    }

#endif

    /** Free some memory **/
    costMatrixExit(ctrl->costmtx);
    eTreeExit(ctrl->etree);

    /** gener the final solverMarix for this processor
     i.e. relative bloc numbering **/
    if(ctrl->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    if(ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "%ld : Genering final SolverMatrix \n", (long)clustnum);
    if (ctrl->iparm[IPARM_DOF_COST] != 0) {
        Dof dofstr;
        dofInit(&dofstr);
        dofConstant(&dofstr, 0, symbmtx->nodenbr, ctrl->iparm[IPARM_DOF_NBR]);
        bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, &dofstr);
        dofExit(&dofstr);
    }
    else{
        bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, dofptr);
    }

    if(ctrl->timer)
    {
        clockStop(timer_current);
        printf("--SolverMatrix computed at time: %g --\n", clockVal(timer_current));
    }

    /*if(ctrl->count_ops)
     {
     printSolverInfo(stderr, solvmtx, dofptr);
     }*/

    if( ctrl->leader == clustnum &&
        ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "** End of Partition & Distribution phase ** \n");

    /** Time end **/
    if(ctrl->timer)
    {
        clockStop(timer_all);
        printf("---- Total execution at time: %g ----\n",clockVal(timer_all));
        set_dparm(ctrl->dparm, DPARM_ANALYZE_TIME, clockVal(timer_all));
    }

    /** Free allocated memory **/
    simuExit(simuctrl, ctrl->clustnbr, ctrl->total_nbcores, ctrl->local_nbctxts);
    eGraphExit(ctrl->egraph);

    if(ctrl->debug)
    {
        setBcofPtr(solvmtx, bcofind);

        if( ctrl->leader == clustnum &&
            ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
            fprintf(stdout, OUT_BLEND_CHKSOLVER);
        solverCheck(solvmtx);
    }

    /***************************************
     * Realloc Memory in a contiguous way  *
     ***************************************/
    printf("Contiguous reallocation of the solverMatrix ...\n");
    solverRealloc(solvmtx, bcofind);
    printf("Done \n");

#ifdef DEBUG_BLEND
    if (leader == clustnum)
        fprintf(stdout, OUT_BLEND_CHKSOLVER);
    if (ctrl->ricar) {
        if (leader == clustnum)
            errorPrintW("No solverMatrix checking in incomplete factorisation.");
    }else {
        solverCheck(solvmtx);
    }
#endif
    if (bcofind != NULL)
        memFree_null(bcofind);
}
