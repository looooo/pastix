/**
 *
 * @file pastix_task_sopalin.c
 *
 *  PaStiX factorization routines
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "isched.h"
#include "spm.h"
#include "bcsc.h"
#include "blend/solver.h"
#include "coeftab.h"
#include "sopalin_data.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_ccores.h"
#include "kernels/pastix_dcores.h"
#include "kernels/pastix_scores.h"

static void (*sopalinFacto[4][4])(pastix_data_t *, sopalin_data_t*) =
{
    { sopalin_spotrf, sopalin_dpotrf, sopalin_cpotrf, sopalin_zpotrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_csytrf, sopalin_zsytrf },
    { sopalin_sgetrf, sopalin_dgetrf, sopalin_cgetrf, sopalin_zgetrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_chetrf, sopalin_zhetrf }
};

static fct_ge2lr_t compressMethod[2][4] =
{
    { core_sge2lr_SVD_interface , core_dge2lr_SVD_interface , core_cge2lr_SVD_interface , core_zge2lr_SVD_interface  },
    { core_sge2lr_RRQR_interface, core_dge2lr_RRQR_interface, core_cge2lr_RRQR_interface, core_zge2lr_RRQR_interface }
};

static fct_rradd_t recompressMethod[2][4] =
{
    { core_srradd_SVD_interface , core_drradd_SVD_interface , core_crradd_SVD_interface , core_zrradd_SVD_interface  },
    { core_srradd_RRQR_interface, core_drradd_RRQR_interface, core_crradd_RRQR_interface, core_zrradd_RRQR_interface }
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Fill the internal block CSC structure.
 *
 * This internal block CSC can be used by the refinement step with or without
 * the preconditioner obtained by the numerical factorization.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal block CSC is filled with entries from
 *          the spm matrix.
 *
 * @param[in] spm
 *          The sparse matrix descriptor that describes problem instance.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_spm2bcsc( pastix_data_t *pastix_data,
                         const pastix_spm_t  *spm )
{
    double time;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_spm2bcsc: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_subtask_spm2bcsc: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_subtask_spm2bcsc: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Compute the norm of A, to scale the epsilon parameter for pivoting
     */
    {
        pastix_data->dparm[ DPARM_A_NORM ] = spmNorm( PastixFrobeniusNorm, spm );
        if (pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO ) {
            pastix_print( 0, 0,
                          "    ||A||_2  =                            %e\n",
                          pastix_data->dparm[ DPARM_A_NORM ] );
        }
    }

    /*
     * Fill in the internal blocked CSC. We consider that if this step is called
     * the spm values have changed so we need to update the blocked csc.
     */
    if (pastix_data->bcsc != NULL)
    {
        bcscExit( pastix_data->bcsc );
        memFree_null( pastix_data->bcsc );
    }

    MALLOC_INTERN( pastix_data->bcsc, 1, pastix_bcsc_t );

    time = bcscInit( spm,
                     pastix_data->ordemesh,
                     pastix_data->solvmatr,
                     ( (pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU)
                       && (! pastix_data->iparm[IPARM_ONLY_REFINE]) ),
                     pastix_data->bcsc );

    if ( pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT ) {
        pastix_print( 0, 0, OUT_BCSC_TIME, time );
    }

    if ( pastix_data->iparm[IPARM_FREE_CSCUSER] ) {
        spmExit( spm );
    }

    /*
     * Invalidate following step, and add current step to the ones performed
     */
    pastix_data->steps &= ~STEP_BCSC2CTAB;
    pastix_data->steps |= STEP_CSC2BCSC;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Fill the internal solver matrix structure.
 *
 * This step is linked with the pastix_subtask_sopalin() since this structure is
 * only used during the numerical factorization. 
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, [IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal solver structure is filled with entries from
 *          the internal block CSC matrix.
 *
 * @param[in] spm
 *          The sparse matrix descriptor that describes problem instance.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_bcsc2ctab( pastix_data_t      *pastix_data,
                          const pastix_spm_t *spm )
{
    Clock timer;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_subtask_bcsc2ctab: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_subtask_bcsc2ctab: All steps from pastix_task_init() to pastix_stask_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    clockStart(timer);

    /* Initialize low-rank parameters */
    pastix_data->solvmatr->lowrank.compress_when       = pastix_data->iparm[IPARM_COMPRESS_WHEN];
    pastix_data->solvmatr->lowrank.compress_method     = pastix_data->iparm[IPARM_COMPRESS_METHOD];
    pastix_data->solvmatr->lowrank.compress_min_width  = pastix_data->iparm[IPARM_COMPRESS_MIN_WIDTH];
    pastix_data->solvmatr->lowrank.compress_min_height = pastix_data->iparm[IPARM_COMPRESS_MIN_HEIGHT];
    pastix_data->solvmatr->lowrank.tolerance           = sqrt(pastix_data->dparm[DPARM_COMPRESS_TOLERANCE]);

    pastix_data->solvmatr->lowrank.core_ge2lr = compressMethod[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][spm->flttype-2];
    pastix_data->solvmatr->lowrank.core_rradd = recompressMethod[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][spm->flttype-2];

    pastix_data->solvmatr->factotype = pastix_data->iparm[IPARM_FACTORIZATION];

    /*
     * Fill in the internal coeftab structure. We consider that if this step is
     * called the bcsc values have changed, or a factorization have already been
     * performed, so we need to update the coeftab arrays.
     */
    if (pastix_data->bcsc != NULL)
    {
        coeftabExit( pastix_data->solvmatr );
    }

    coeftabInit( pastix_data,
                 spm->flttype == PastixPattern,
                 pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU );

#if defined(PASTIX_WITH_PARSEC)
    {
        sparse_matrix_desc_t *sdesc = pastix_data->solvmatr->parsec_desc;
        int mtxtype = ( pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU ) ? PastixGeneral : PastixHermitian;
        if ( sdesc != NULL ) {
            sparse_matrix_destroy( sdesc );
        }
        else {
            sdesc = (sparse_matrix_desc_t*)malloc(sizeof(sparse_matrix_desc_t));
        }

        /* Create the matrix descriptor */
        sparse_matrix_init( sdesc, pastix_data->solvmatr,
                            pastix_size_of( spm->flttype ), mtxtype,
                            1, 0 );
        pastix_data->solvmatr->parsec_desc = sdesc;
    }
#endif

    clockStop(timer);
    if (pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
        pastix_print( 0, 0, OUT_COEFTAB_TIME,
                      clockVal(timer) );
    }

    /* Invalidate following step, and add current step to the ones performed */
    pastix_data->steps &= ~STEP_NUMFACT;
    pastix_data->steps |= STEP_BCSC2CTAB;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_numfact
 *
 * @brief Factorize the given problem using Cholesky or LU decomposition.
 *
 * The user can call the pastix_task_solve() to obtain the solution.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_COMPRESS_WHEN.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the solver matrix structure stores the factorization of the
 *          given problem.
 *
 * @param[in] spm
 *          The sparse matrix descriptor that describes problem instance.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_subtask_sopalin( pastix_data_t *pastix_data,
                        const pastix_spm_t  *spm )
{
    sopalin_data_t  sopalin_data;
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
    if (spm == NULL) {
        errorPrint("pastix_task_sopalin: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm   = pastix_data->iparm;
    procnum = pastix_data->inter_node_procnum;

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
        pastix_print(procnum, 0, OUT_STEP_SOPALIN,
                     pastixFactotypeStr( iparm[IPARM_FACTORIZATION] ) );
    }

    /* Prepare the sopalin_data structure */
    {
        sopalin_data.solvmtx = pastix_data->solvmatr;

        /* TODO: might change the behavior: if the user wants a ratio of the norm, it could compute it himself */
        if ( pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] < 0. ) {
            sopalin_data.diagthreshold = - pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ];
        }
        else if ( pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] == 0. ) {
            if ( (spm->flttype == PastixFloat) || (spm->flttype == PastixComplex32) )
                sopalin_data.diagthreshold = 1e-7  * pastix_data->dparm[DPARM_A_NORM];
            else
                sopalin_data.diagthreshold = 1e-15 * pastix_data->dparm[DPARM_A_NORM];
        }
        else {
            sopalin_data.diagthreshold = pastix_data->dparm[ DPARM_EPSILON_MAGN_CTRL ] * pastix_data->dparm[DPARM_A_NORM];
        }
    }

    sbackup = solverBackupInit( pastix_data->solvmatr );
    pastix_data->solvmatr->restore = 2;
    {
        void (*factofct)( pastix_data_t *, sopalin_data_t *);
        double timer, flops;

        factofct = sopalinFacto[ pastix_data->iparm[IPARM_FACTORIZATION] ][spm->flttype-2];
        assert(factofct);

        clockStart(timer);
        factofct( pastix_data, &sopalin_data );
        clockStop(timer);

        flops = pastix_data->dparm[DPARM_FACT_FLOPS] / clockVal(timer);
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
            pastix_print( 0, 0, OUT_SOPALIN_TIME,
                          clockVal(timer),
                          printflopsv( flops ), printflopsu( flops ) );
        }
    }
    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbol.eps", "w");
        solverDraw(pastix_data->solvmatr,
                   stream,
                   iparm[IPARM_VERBOSE]);
        fclose(stream);
    }
#endif

    if ( (pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) &&
         (pastix_data->iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever) )
    {
        /* Compute the memory gain */
        coeftabMemory[spm->flttype-2]( pastix_data->solvmatr );
    }

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_BCSC2CTAB |
                             STEP_SOLVE     |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    return EXIT_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Perform all the numerical factorization steps:
 * fill the internal block CSC and the solver matrix structures,
 * then apply the factorization step.
 *
 * The user can call the pastix_task_solve() to obtain the solution.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal block CSC is filled with entries from the
 *          spm matrix and the solver matrix structure stores the factorization
 *          of the given problem.
 *
 * @param[in] spm
 *          The sparse matrix descriptor that describes problem instance.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastix_task_numfact( pastix_data_t *pastix_data,
                     const pastix_spm_t  *spm )
{
    /* #ifdef PASTIX_WITH_MPI */
    /*     MPI_Comm       pastix_comm = pastix_data->inter_node_comm; */
    /* #endif */
    pastix_int_t  procnum;
    pastix_int_t *iparm;
    /*     double        *dparm    = pastix_data->dparm; */
    /*     SolverMatrix  *solvmatr = pastix_data->solvmatr; */
    int rc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_sopalin: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_task_sopalin: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_ANALYSE) ) {
        errorPrint("pastix_task_sopalin: All steps from pastix_task_init() to pastix_task_blend() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm   = pastix_data->iparm;
    procnum = pastix_data->inter_node_procnum;

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
        pastix_print(procnum, 0, OUT_STEP_SOPALIN,
                     pastixFactotypeStr( iparm[IPARM_FACTORIZATION] ) );
    }

    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        rc = pastix_subtask_spm2bcsc( pastix_data, spm );
        if (rc != PASTIX_SUCCESS)
            return rc;
    }

    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        rc = pastix_subtask_bcsc2ctab( pastix_data, spm );
        if (rc != PASTIX_SUCCESS)
            return rc;
    }

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        rc = pastix_subtask_sopalin( pastix_data, spm );
        if (rc != PASTIX_SUCCESS)
            return rc;
    }

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_BCSC2CTAB |
                             STEP_SOLVE     |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    return EXIT_SUCCESS;
}
