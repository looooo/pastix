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
#include "param_blend.h"
#include "blendctrl.h"
#include "splitpart.h"
#include "write_ps.h"
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
 *   clustnbr   - Number of MPI processes.
 *   thrdlocnbr - Number of threads.
 *   cudanbr    - Number of cuda devices.
 *   clustnum   - Processor ID number.
 *   option     - Blend parameters.
 *   dofptr     -
 */
void solverBlend(SolverMatrix *solvmtx,
                 SymbolMatrix *symbmtx,
                 int           clustnbr,
                 int           thrdlocnbr,
                 int           cudanbr,
                 int           clustnum,
                 BlendParam   *option,
                 const Dof    *dofptr)
{
    BlendCtrl *ctrl;
    SimuCtrl  *simuctrl;
    double     timer_all     = 0.;
    double     timer_current = 0.;
    pastix_int_t       *bcofind       = NULL;

    /* Initialisation of the control structure */
    MALLOC_INTERN(ctrl, 1, BlendCtrl);
    blendCtrlInit(ctrl, clustnbr, thrdlocnbr, cudanbr, clustnum, option);

    /* Check parameters */
    if(clustnum >= clustnbr)
    {
        errorPrint("solverBlend parameter clustnum(%ld) is greater"
                   " than clustnbr(%ld).",
                   (long) clustnum, (long) clustnbr);
        EXIT(MOD_BLEND,INTERNAL_ERR);
    }
    if(ctrl->option->ooc)
    {
        /* OOC works only with 1D structures */
        if ((ctrl->option->leader == clustnum) &&
            (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
            fprintf(stdout, "Force 1D distribution because of OOC \n");
        ctrl->option->ratiolimit = INTVALMAX;
    }

    clockStart(timer_all);

    pastix_print( clustnum, 0,
                  OUT_CLUSTNBR "" OUT_PROCNBR "" OUT_THRDNBR,
                  (long)clustnbr, (long)ctrl->proclocnbr, (long)ctrl->thrdlocnbr);

    /* Get the symbmtx ready */
    {
        symbolRealloc(symbmtx);
        symbolBase(symbmtx, 0);
        symbolCheck(symbmtx);

        /* Rustine */
#define RUSTINE
#ifdef RUSTINE
        symbolRustine(symbmtx, symbmtx);
#endif
    }

    /** Verify the coherence of the initial symbol matrix **/
    if(ctrl->option->debug)
    {
        pastix_print( clustnum, 0, OUT_BLEND_CHKSMBMTX );
        symbolCheck(symbmtx);
    }

    if(ctrl->option->count_ops && ctrl->option->leader == clustnum)
        symbCost(option->iparm, option->dparm, symbmtx, dofptr);

    /* build the elimination graph from the symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_ELIMGRAPH );
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);
        eGraphInit(ctrl->egraph);
        eGraphBuild(ctrl->egraph, symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Graph build at time: %g --\n", clockVal(timer_current) );
    }

    /* Build the elimination tree from the symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_ELIMTREE );
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->etree, 1, EliminTree);
        eTreeInit(ctrl->etree);
        eTreeBuild(ctrl->etree, symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Tree build at time: %g --\n", clockVal(timer_current));
    }

    /* Build the cost matrix from the symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_COSTMATRIX );
        clockStart(timer_current);

        MALLOC_INTERN(ctrl->costmtx, 1, CostMatrix);
        costInit(ctrl->costmtx);
        costMatrixBuild(ctrl->costmtx, symbmtx, dofptr);

        /* Compute costs for all nodes */
        subtreeUpdateCost(eTreeRoot(ctrl->etree), ctrl->costmtx, ctrl->etree);
        if(ctrl->option->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            fprintf(stdout, "Total cost of the elimination tree %g \n",
                    ctrl->costmtx->cblktab[eTreeRoot(ctrl->etree)].subtree);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Cost Matrix build at time: %g --\n", clockVal(timer_current));
    }

    /*
     * Partitioning of the initial symbolic factorization and processing of
     * candidate processors group for each colum bloc
     */
    {
        pastix_print( clustnum, 0, "Spliting initial partition \n" );
        clockStart(timer_current);

        splitPart(symbmtx, ctrl, dofptr);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Split build at time: %g --\n", clockVal(timer_current));
    }

#if defined(PASTIX_DEBUG_BLEND)
    /** Verify the coherence of the new symbol matrix **/
    symbolCheck(symbmtx);
#endif

    //TODO
    /* if ( (ctrl->option->leader == clustnum) && (ctrl->option->tracegen == 1)) */
    /*   { */
    /*     FILE *out; */
    /*     OUT_OPENFILEINDIR(ctrl->option->iparm, out, "elimintree.dot", "w"); */
    /*     eTreeGenDot(ctrl->etree, out); */
    /*     OUT_CLOSEFILEINDIR(out); */
    /*   } */

    if(ctrl->option->count_ops && (ctrl->option->leader == clustnum))
        symbCost(option->iparm, option->dparm, symbmtx, dofptr);

    pastix_print( clustnum, 0, "** New Partition: cblknbr=  %ld     bloknbr=  %ld     ratio=%f ** \n",
                  (long)symbmtx->cblknbr, (long)symbmtx->bloknbr,
                  (float)symbmtx->bloknbr/(float)symbmtx->cblknbr);

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbolblend.eps", "w");
        symbolDraw(symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    /* Re-build the new elimination graph from the new symbolic partition */
    {
        pastix_print( clustnum, 0, OUT_BLEND_ELIMGRAPH2 );
        clockStart(timer_current);

        eGraphExit(ctrl->egraph);
        MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);
        eGraphInit(ctrl->egraph);
        eGraphBuild(ctrl->egraph, symbmtx);

        clockStop(timer_current);
        pastix_print( clustnum, 0, "--Graph build at time: %g --\n", clockVal(timer_current) );
    }

    if(ctrl->option->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    /* initialize simu structure control */
    MALLOC_INTERN(simuctrl, 1, SimuCtrl);
    simuInit(simuctrl, symbmtx, ctrl->clustnbr, ctrl->procnbr,
             symbmtx->cblknbr, symbmtx->bloknbr, ctrl->candtab);

    /* Build tasks */
    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_TASKGRAPH);
    taskBuild(simuctrl, symbmtx, ctrl->candtab, dofptr, ctrl);

    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_NBTASK, (long)simuctrl->tasknbr);
    if(ctrl->option->timer)
    {
        clockStop(timer_current);
        printf("--Task built at time: %g --\n", clockVal(timer_current));
    }

    /* Distribution Phase */
    if(ctrl->option->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

#ifdef DEBUG_BLEND
    ASSERT(check_candidat(symbmtx, ctrl)>=0,MOD_BLEND);
#endif

    if((ctrl->option->leader == clustnum) &&
       (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
        fprintf(stdout, OUT_BLEND_DISTPART);
    distribPart(symbmtx, simuctrl, ctrl, dofptr);

    if(ctrl->option->timer)
    {
        clockStop(timer_current);
        printf("--Distribution computed at time: %g --\n", clockVal(timer_current));
    }

#ifdef PASTIX_DYNSCHED /* 2 eme passe de splitpart */

    if(ctrl->option->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    /** repartitioning of the initial symbolic factorization
     and processing of candidate processors group for
     each colum bloc **/

    splitPartLocal(ctrl, simuctrl, symbmtx, dofptr);

    if(ctrl->option->timer)
    {
        clockStop(timer_current);
        printf("--Split build at time: %g --\n", clockVal(timer_current));
    }

#endif

    /** Free some memory **/
    costExit(ctrl->costmtx);
    eTreeExit(ctrl->etree);

    /** gener the final solverMarix for this processor
     i.e. relative bloc numbering **/
    if(ctrl->option->timer)
    {
        clockInit(timer_current);
        clockStart(timer_current);
    }

    if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "%ld : Genering final SolverMatrix \n", (long)clustnum);
    if (ctrl->option->iparm[IPARM_DOF_COST] != 0) {
        Dof dofstr;
        dofInit(&dofstr);
        dofConstant(&dofstr, 0, symbmtx->nodenbr, ctrl->option->iparm[IPARM_DOF_NBR]);
        bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, &dofstr);
        dofExit(&dofstr);
    }
    else{
        bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, dofptr);
    }

    if(ctrl->option->timer)
    {
        clockStop(timer_current);
        printf("--SolverMatrix computed at time: %g --\n", clockVal(timer_current));
    }

    /*if(ctrl->option->count_ops)
     {
     printSolverInfo(stderr, solvmtx, dofptr);
     }*/

    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "** End of Partition & Distribution phase ** \n");

    /** Time end **/
    if(ctrl->option->timer)
    {
        clockStop(timer_all);
        printf("---- Total execution at time: %g ----\n",clockVal(timer_all));
        set_dparm(option->dparm, DPARM_ANALYZE_TIME, clockVal(timer_all));
    }

    /** Free allocated memory **/
    simuExit(simuctrl, ctrl->clustnbr, ctrl->procnbr, ctrl->bublnbr);
    eGraphExit(ctrl->egraph);

    if(ctrl->option->debug)
    {
        setBcofPtr(solvmtx, bcofind);

        if( ctrl->option->leader == clustnum &&
            ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
            fprintf(stdout, OUT_BLEND_CHKSOLVER);
        solverCheck(solvmtx);
    }

    /***************************************
     * Realloc Memory in a contiguous way  *
     ***************************************/
    blendCtrlExit(ctrl);
    printf("Contiguous reallocation of the solverMatrix ...\n");
    solverRealloc(solvmtx, bcofind);
    printf("Done \n");

#ifdef DEBUG_BLEND
    if (leader == clustnum)
        fprintf(stdout, OUT_BLEND_CHKSOLVER);
    if (option->ricar) {
        if (leader == clustnum)
            errorPrintW("No solverMatrix checking in incomplete factorisation.");
    }else {
        solverCheck(solvmtx);
    }
#endif
    if (bcofind != NULL)
        memFree_null(bcofind);
}
