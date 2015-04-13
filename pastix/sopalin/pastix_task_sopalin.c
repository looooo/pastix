/**
 *
 * @file pastix_task_sopalin.c
 *
 *  PaStiX factorization routines
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
#include "csc.h"
#include "bcsc.h"
#include "sopalin_data.h"

void
coeftabInit( const SolverMatrix  *datacode,
             const pastix_bcsc_t *bcsc,
             pastix_int_t         fakefillin,
             pastix_int_t         factoLU );

/**
 *******************************************************************************
 *
 * @ingroup pastix_sopalin
 * @ingroup pastix
 *
 * pastix_task_sopalin - Factorize the given problem using Cholesky(-Crout) or
 * LU decomposition.
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
 * @param[in,out] csc
 *          ...
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
pastix_task_sopalin( pastix_data_t *pastix_data,
                     pastix_csc_t  *csc )
{
    SolverBackup_t *sbackup;
/* #ifdef PASTIX_WITH_MPI */
/*     MPI_Comm       pastix_comm = pastix_data->inter_node_comm; */
/* #endif */
    pastix_int_t  procnum;
    pastix_int_t *iparm;
/*     double        *dparm    = pastix_data->dparm; */
/*     SolverMatrix  *solvmatr = pastix_data->solvmatr; */

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_sopalin: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (csc == NULL) {
        errorPrint("pastix_task_sopalin: wrong csc parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm   = pastix_data->iparm;
    procnum = pastix_data->inter_node_procnum;

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            pastix_print(procnum, 0, "%s", OUT_STEP_NUMFACT_LU);
            break;
        case API_FACT_LLT:
            pastix_print(procnum, 0, "%s", OUT_STEP_NUMFACT_LLT);
            break;
        case API_FACT_LDLH:
            pastix_print(procnum, 0, "%s", OUT_STEP_NUMFACT_LDLH);
            break;
        case API_FACT_LDLT:
        default:
            pastix_print(procnum, 0, "%s", OUT_STEP_NUMFACT_LDLT);
        }
    }

    /**
     * Fill in the internal blocked CSC. We consider that if this step is called
     * the csc values have changed so we need to update the blocked csc.
     */
    if (pastix_data->bcsc != NULL)
    {
        bcscExit( pastix_data->bcsc );
        memFree_null( pastix_data->bcsc );
    }

    MALLOC_INTERN( pastix_data->bcsc, 1, pastix_bcsc_t );

    /**
     * Check that IPARM_ONLY_RAFF is not set, because we don't need to
     * initialize the A^t array in this case and we should have never arrived
     * here.
     */
    assert( !iparm[IPARM_ONLY_RAFF] );

    bcscInit( csc,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU,
              pastix_data->bcsc );

    if ( iparm[IPARM_FREE_CSCUSER] ) {
        spmExit( csc );
    }

    coeftabInit( pastix_data->solvmatr,
                 pastix_data->bcsc,
                 csc->flttype == PastixPattern,
                 csc->mtxtype == PastixGeneral );
    sbackup = solverBackupInit( pastix_data->solvmatr );

    {
        sopalin_data_t sopalin_data;
        double timer;

        clockStart(timer);
        sopalin_data.solvmtx = pastix_data->solvmatr;

        pastix_static_dpotrf( &sopalin_data );
        clockStop(timer);
        pastix_print( 0, 0, "--Factorization time: %g --\n", clockVal(timer));

        //pastix_static_dsytrf( &sopalin_data );
        //coeftab_ddump( pastix_data->solvmatr, "AfterFacto" );
    }

        /* sopalinInit( pastix_data, sopalin_data, 1 ); */
        /* sopalinFacto( pastix_data, sopalin_data ); */
        /* sopalinExit( sopalin_data ); */

    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_SOLVE  |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    iparm[IPARM_START_TASK]++;

    return EXIT_SUCCESS;
}
