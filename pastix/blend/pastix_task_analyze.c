/**
 *
 * @file pastix_task_analyze.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "order.h"
#include "perf.h"
#include "elimin.h"
#include "cost.h"
#include "cand.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "blend.h"
#include "dof.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_analyze
 * @ingroup pastix
 *
 * pastix_task_blend - Computes the structural information required to factorize
 * and solve the problem.
 *
 * ...
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, ...
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, ...
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
// TODO: need some cleanup arount the printf, and solverBlend could be
// integrated directly into this function.
int
pastix_task_blend(pastix_data_t *pastix_data)
{
    Dof            dofstr;
    BlendCtrl      ctrl;
#ifdef PASTIX_WITH_MPI
    MPI_Comm       pastix_comm = pastix_data->inter_node_comm;
#endif
    pastix_int_t   procnbr  = pastix_data->inter_node_procnbr;
    pastix_int_t   procnum  = pastix_data->inter_node_procnum;
    pastix_int_t  *iparm    = pastix_data->iparm;
    double        *dparm    = pastix_data->dparm;
    SolverMatrix  *solvmatr;

    if ( !(pastix_data->steps & STEP_SYMBFACT) ) {
        errorPrint("pastix_task_symbfact: pastix_task_init() has to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print( procnum, 0, "%s", OUT_STEP_BLEND );

    if ( pastix_data->solvmatr != NULL ) {
        solverExit( pastix_data->solvmatr );
        memFree_null( pastix_data->solvmatr );
    }
    solvmatr = (SolverMatrix*)malloc(sizeof(SolverMatrix));
    pastix_data->solvmatr = solvmatr;

    dofInit(&dofstr);
    dofConstant(&dofstr, 0, pastix_data->symbmtx->nodenbr,
                ((iparm[IPARM_DOF_COST] == 0)?iparm[IPARM_DOF_NBR]:iparm[IPARM_DOF_COST]));

    blendCtrlInit( &ctrl, procnum, procnbr,
                   iparm[IPARM_THREAD_NBR],
                   iparm[IPARM_THREAD_NBR],
                   iparm, dparm );

#ifdef FORCE_NOSMP
    iparm[IPARM_THREAD_NBR] = 1;
#endif

    solverBlend( &ctrl, solvmatr, pastix_data->symbmtx, &dofstr );
    blendCtrlExit(&ctrl);

    symbolPrintStats(pastix_data->symbmtx);

    symbolExit(pastix_data->symbmtx);
    memFree_null(pastix_data->symbmtx);

    if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
    {
        iparm[IPARM_NNZEROS]       *=2;
        dparm[DPARM_PRED_FACT_TIME]*=2.;
    }
    dparm[DPARM_SOLV_FLOPS] = (double)iparm[IPARM_NNZEROS]; /* number of operations for solve */

    iparm[IPARM_NNZEROS_BLOCK_LOCAL] = solvmatr->coefnbr;

    /* Affichage */
    dparm[DPARM_FILL_IN] = dparm[DPARM_FILL_IN]
        * (double)(iparm[IPARM_NNZEROS]/(iparm[IPARM_DOF_NBR]*iparm[IPARM_DOF_NBR]));

    if ((procnum==0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
    {
        fprintf(stdout, TIME_TO_ANALYSE,       (double)dparm[DPARM_ANALYZE_TIME]);
        fprintf(stdout, NNZERO_WITH_FILLIN_TH, (long)iparm[IPARM_NNZEROS]);
        fprintf(stdout, OUT_FILLIN_TH,         (double)dparm[DPARM_FILL_IN]);
        if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
            fprintf(stdout,NUMBER_OP_LU,     (double)dparm[DPARM_FACT_FLOPS]);
        else
            fprintf(stdout, NUMBER_OP_LLT,    (double)dparm[DPARM_FACT_FLOPS]);
        fprintf(stdout,TIME_FACT_PRED,PERF_MODEL,     (double)dparm[DPARM_PRED_FACT_TIME]);

    }
    if ((iparm[IPARM_VERBOSE] > API_VERBOSE_NO))
    {
        pastix_int_t solversize = sizeofsolver(solvmatr, iparm);

        fprintf(stdout,SOLVMTX_WITHOUT_CO,  (int)procnum, (double)MEMORY_WRITE(solversize),
                MEMORY_UNIT_WRITE(solversize));
        fprintf(stdout, NNZERO_WITH_FILLIN, (int)procnum, (long)iparm[IPARM_NNZEROS_BLOCK_LOCAL]);
    }
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    {
        pastix_int_t sizeL = solvmatr->coefnbr;
        pastix_int_t sizeG = 0;

        MPI_Reduce(&sizeL, &sizeG, 1, PASTIX_MPI_INT, MPI_MAX, 0, pastix_comm);

        if (procnum == 0)
        {
            sizeG *= sizeof(pastix_complex64_t);
            if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
                sizeG *= 2;

            fprintf(stdout, OUT_COEFSIZE, (double)MEMORY_WRITE(sizeG),
                    MEMORY_UNIT_WRITE(sizeG));
        }
    }

    if (1)
    {
        FILE *file = fopen("solvergen", "w");
        solverSave( solvmatr, file );
        fclose(file);
    }

    /* Invalidate following steps, and add analyse step to the ones performed */
    pastix_data->steps &= ~( STEP_NUMFACT |
                             STEP_SOLVE   |
                             STEP_REFINE  );
    pastix_data->steps |= STEP_ANALYSE;

    iparm[IPARM_START_TASK]++;

    return PASTIX_SUCCESS;
}
