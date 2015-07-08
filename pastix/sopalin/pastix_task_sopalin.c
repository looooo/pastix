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
#include "isched.h"
#include "csc.h"
#include "bcsc.h"
#include "sopalin_data.h"

static void (*sopalinFacto[4][4])(sopalin_data_t*) =
{
    { sopalin_spotrf, sopalin_dpotrf, sopalin_cpotrf, sopalin_zpotrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_csytrf, sopalin_zsytrf },
    { sopalin_sgetrf, sopalin_dgetrf, sopalin_cgetrf, sopalin_zgetrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_chetrf, sopalin_zhetrf }
};

void
coeftabInit( const SolverMatrix  *datacode,
             const pastix_bcsc_t *bcsc,
             pastix_int_t         fakefillin,
             pastix_int_t         factoLU );

int
pastix_subtask_csc2bcsc( pastix_data_t *pastix_data,
                         pastix_csc_t  *csc )
{
    /**
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_csc2bcsc: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (csc == NULL) {
        errorPrint("pastix_subtask_csc2bcsc: wrong csc parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_subtask_csc2bcsc: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
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

    bcscInit( csc,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              ( (pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU)
                && (! pastix_data->iparm[IPARM_ONLY_RAFF]) ),
              pastix_data->bcsc );

    if ( pastix_data->iparm[IPARM_FREE_CSCUSER] ) {
        spmExit( csc );
    }

    /* Invalidate following step, and add current step to the ones performed */
    pastix_data->steps &= ~STEP_BCSC2CTAB;
    pastix_data->steps |= STEP_CSC2BCSC;

    return PASTIX_SUCCESS;
}


int
pastix_subtask_bcsc2ctab( pastix_data_t *pastix_data,
                          pastix_csc_t  *csc )
{
    /**
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (csc == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong csc parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_subtask_bcsc2ctab: All steps from pastix_task_init() to pastix_stask_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    coeftabInit( pastix_data->solvmatr,
                 pastix_data->bcsc,
                 csc->flttype == PastixPattern,
                 pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU );

    /* Invalidate following step, and add current step to the ones performed */
    pastix_data->steps &= ~STEP_NUMFACT;
    pastix_data->steps |= STEP_BCSC2CTAB;

    return PASTIX_SUCCESS;
}


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
    extern isched_t *scheduler;
    sopalin_data_t  sopalin_data;
    SolverBackup_t *sbackup;
/* #ifdef PASTIX_WITH_MPI */
/*     MPI_Comm       pastix_comm = pastix_data->inter_node_comm; */
/* #endif */
    pastix_int_t  procnum;
    pastix_int_t *iparm;
    int rc;
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


    /* Compute the norm of A, to scale the epsilon parameter for pivoting */
    {
        pastix_print( 0, 0, "-- ||A||_2  =                                   " );
        pastix_data->dparm[ DPARM_A_NORM ] = spmNorm( PastixFrobeniusNorm, csc );
        pastix_print( 0, 0, "%lg\n", pastix_data->dparm[ DPARM_A_NORM ] );
    }

    rc = pastix_subtask_csc2bcsc( pastix_data, csc );
    if (rc != PASTIX_SUCCESS)
        return rc;

    rc = pastix_subtask_bcsc2ctab( pastix_data, csc );
    if (rc != PASTIX_SUCCESS)
        return rc;

    /* Prepare the sopalin_data structure */
    {
        sopalin_data.solvmtx = pastix_data->solvmatr;

        /* TODO: might change the behavior: if the user wants a ratio of the norm, it could compute it himself */
        if ( pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] < 0. ) {
            sopalin_data.diagthreshold = - pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ];
        }
        else if ( pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] == 0. ) {
            if ( (csc->flttype == PastixFloat) || (csc->flttype == PastixComplex32) )
                sopalin_data.diagthreshold = 1e-7  * pastix_data->dparm[DPARM_A_NORM];
            else
                sopalin_data.diagthreshold = 1e-15 * pastix_data->dparm[DPARM_A_NORM];
        }
        else {
            sopalin_data.diagthreshold = pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] * pastix_data->dparm[DPARM_A_NORM];
        }
    }
    sopalin_data.sched = scheduler;

    sbackup = solverBackupInit( pastix_data->solvmatr );
    pastix_data->solvmatr->restore = 2;
    {
        void (*factofct)( sopalin_data_t *);
        double timer;

        factofct = sopalinFacto[ pastix_data->iparm[IPARM_FACTORIZATION] ][csc->flttype-2];
        assert(sopalinFacto);

        clockStart(timer);
        factofct( &sopalin_data );
        clockStop(timer);
        pastix_print( 0, 0, OUT_TIME_FACT, clockVal(timer) );

        //pastix_static_dsytrf( &sopalin_data );
        //coeftab_ddump( pastix_data->solvmatr, "AfterFacto" );
        /* sopalinInit( pastix_data, sopalin_data, 1 ); */
        /* sopalinFacto( pastix_data, sopalin_data ); */
        /* sopalinExit( sopalin_data ); */
    }
    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_SOLVE  |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    iparm[IPARM_START_TASK]++;

    return EXIT_SUCCESS;
}
