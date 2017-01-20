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
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_analyze
 * @ingroup pastix
 *
 * pastix_task_blend - Computes the structural information required to factorize
 * and solve the problem. It requires an ordering structure, as well as the
 * symbolic factorization structure. It results a solver structure that contains
 * all informations dependent of the architecture and problem to solve the
 * problem efficiently.
 * On exit, the symbol structure is destroyed and only local uncompressed
 * information is stored in the solver structure.
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
int
pastix_task_blend(pastix_data_t *pastix_data)
{
    BlendCtrl      ctrl;
    pastix_int_t   procnbr, procnum;
    pastix_int_t  *iparm;
    double        *dparm;
    Order         *ordeptr;
    SymbolMatrix  *symbptr;
    SolverMatrix  *solvptr;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_analyze: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_SYMBFACT) ) {
        errorPrint("pastix_task_analyze: pastix_task_symbfact() has to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }


    /*
     * Free graph structure
     */
    if (pastix_data->graph != NULL) {
        graphExit( pastix_data->graph );
        memFree_null( pastix_data->graph );
    }


    iparm   = pastix_data->iparm;
    dparm   = pastix_data->dparm;
    procnbr = pastix_data->inter_node_procnbr;
    procnum = pastix_data->inter_node_procnum;
    ordeptr = pastix_data->ordemesh;
    symbptr = pastix_data->symbmtx;

    if (ordeptr == NULL) {
        errorPrint("pastix_task_analyze: the pastix_data->ordemesh field has not been initialized, pastix_task_order should be called first");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (symbptr == NULL) {
        errorPrint("pastix_task_analyze: the pastix_data->symbmtx has not been initialized, pastix_task_symbfact should be called first");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (symbptr->dof < 1) {
        errorPrint("pastix_task_analyze: Dof number has not been correctly initialized");
        return PASTIX_ERR_BADPARAMETER;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print( procnum, 0, OUT_STEP_BLEND );

    /* Cleanup the solver structure if we already computed it */
    if ( pastix_data->solvmatr != NULL ) {
        solverExit( pastix_data->solvmatr );
        memFree_null( pastix_data->solvmatr );
    }
    solvptr = (SolverMatrix*)malloc(sizeof(SolverMatrix));
    pastix_data->solvmatr = solvptr;

    blendCtrlInit( &ctrl, procnum, procnbr,
                   iparm, dparm );

    solverBlend( &ctrl, solvptr, symbptr );

    blendCtrlExit(&ctrl);

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        symbolPrintStats( pastix_data->symbmtx );

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbol.eps", "w");
        symbolDraw(pastix_data->symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    /* Symbol is not used anymore */
    symbolExit(pastix_data->symbmtx);
    memFree_null(pastix_data->symbmtx);

    /* Computes and print statistics */
    {
#ifdef PASTIX_WITH_MPI
        MPI_Comm pastix_comm = pastix_data->inter_node_comm;
#endif
        if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
        {
            iparm[IPARM_NNZEROS]        *= 2;
            dparm[DPARM_PRED_FACT_TIME] *= 2.;
        }
        dparm[DPARM_SOLV_FLOPS] = (double)iparm[IPARM_NNZEROS]; /* number of operations for solve */

        iparm[IPARM_NNZEROS_BLOCK_LOCAL] = solvptr->coefnbr;

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
            fprintf(stdout, NNZERO_WITH_FILLIN, (int)procnum, (long)iparm[IPARM_NNZEROS_BLOCK_LOCAL]);
        }
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        {
            pastix_int_t sizeL = solvptr->coefnbr;
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
    }

    /* Backup the solver for debug */
    if (0)
    {
        FILE *file = fopen("solvergen", "w");
        solverSave( solvptr, file );
        fclose(file);
    }

    /* Invalidate following steps, and add analyze step to the ones performed */
    pastix_data->steps &= ~( STEP_NUMFACT |
                             STEP_SOLVE   |
                             STEP_REFINE  );
    pastix_data->steps |= STEP_ANALYSE;

    iparm[IPARM_START_TASK]++;

    return PASTIX_SUCCESS;
}
