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
#include "spm.h"
#include "bcsc.h"
#include "blend/solver.h"
#include "coeftab.h"
#include "sopalin_data.h"

void core_sge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_dge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_cge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_zge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_sge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_dge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_cge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_zge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );

int core_srradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_drradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_crradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_zrradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_srradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_drradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_crradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_zrradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);

static void (*sopalinFacto[4][4])(pastix_data_t *, sopalin_data_t*) =
{
    { sopalin_spotrf, sopalin_dpotrf, sopalin_cpotrf, sopalin_zpotrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_csytrf, sopalin_zsytrf },
    { sopalin_sgetrf, sopalin_dgetrf, sopalin_cgetrf, sopalin_zgetrf },
    { sopalin_ssytrf, sopalin_dsytrf, sopalin_chetrf, sopalin_zhetrf }
};

static void (*compressMethod[2][4])(double , pastix_int_t , pastix_int_t ,
                                    void *, pastix_int_t ,
                                    void * ) =
{
    { &core_sge2lr_SVD_interface , &core_dge2lr_SVD_interface , &core_cge2lr_SVD_interface , &core_zge2lr_SVD_interface  },
    { &core_sge2lr_RRQR_interface, &core_dge2lr_RRQR_interface, &core_cge2lr_RRQR_interface, &core_zge2lr_RRQR_interface }
};

static int (*recompressMethod[2][4])(double , int, void *,
                                      pastix_int_t, pastix_int_t, const pastix_lrblock_t *,
                                      pastix_int_t, pastix_int_t,       pastix_lrblock_t *,
                                      pastix_int_t, pastix_int_t ) =
{
    { &core_srradd_SVD_interface , &core_drradd_SVD_interface , &core_crradd_SVD_interface , &core_zrradd_SVD_interface  },
    { &core_srradd_RRQR_interface, &core_drradd_RRQR_interface, &core_crradd_RRQR_interface, &core_zrradd_RRQR_interface }
};

int
pastix_subtask_spm2bcsc( pastix_data_t *pastix_data,
                         pastix_spm_t  *spm )
{
    /**
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

    /**
     * Compute the norm of A, to scale the epsilon parameter for pivoting
     */
    {
        pastix_print( 0, 0, "-- ||A||_2  =                                   " );
        pastix_data->dparm[ DPARM_A_NORM ] = spmNorm( PastixFrobeniusNorm, spm );
        pastix_print( 0, 0, "%lg\n", pastix_data->dparm[ DPARM_A_NORM ] );
    }

    /**
     * Fill in the internal blocked CSC. We consider that if this step is called
     * the spm values have changed so we need to update the blocked csc.
     */
    if (pastix_data->bcsc != NULL)
    {
        bcscExit( pastix_data->bcsc );
        memFree_null( pastix_data->bcsc );
    }

    MALLOC_INTERN( pastix_data->bcsc, 1, pastix_bcsc_t );

    bcscInit( spm,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              ( (pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU)
                && (! pastix_data->iparm[IPARM_ONLY_RAFF]) ),
              pastix_data->bcsc );

    if ( pastix_data->iparm[IPARM_FREE_CSCUSER] ) {
        spmExit( spm );
    }

    /**
     * Invalidate following step, and add current step to the ones performed
     */
    pastix_data->steps &= ~STEP_BCSC2CTAB;
    pastix_data->steps |= STEP_CSC2BCSC;

    return PASTIX_SUCCESS;
}


int
pastix_subtask_bcsc2ctab( pastix_data_t *pastix_data,
                          pastix_spm_t  *spm )
{
    /**
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

    /* Copy the compress_size parameter into the SolverMatrix structure */
    pastix_data->solvmatr->lowrank.compress_when   = pastix_data->iparm[IPARM_COMPRESS_WHEN];
    pastix_data->solvmatr->lowrank.compress_method = pastix_data->iparm[IPARM_COMPRESS_METHOD];
    pastix_data->solvmatr->lowrank.compress_size   = pastix_data->iparm[IPARM_COMPRESS_SIZE];
    pastix_data->solvmatr->lowrank.tolerance       = pastix_data->dparm[DPARM_COMPRESS_TOLERANCE];

    pastix_data->solvmatr->lowrank.core_ge2lr = compressMethod[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][spm->flttype-2];
    pastix_data->solvmatr->lowrank.core_rradd = recompressMethod[ pastix_data->iparm[IPARM_COMPRESS_METHOD] ][spm->flttype-2];

    /**
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
 * @param[in,out] spm
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
                     pastix_spm_t  *spm )
{
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
        double timer;

        factofct = sopalinFacto[ pastix_data->iparm[IPARM_FACTORIZATION] ][spm->flttype-2];
        assert(factofct);

        clockStart(timer);
        factofct( pastix_data, &sopalin_data );
        clockStop(timer);
        pastix_print( 0, 0, OUT_TIME_FACT, clockVal(timer) );
    }
    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    FILE *stream;
    PASTIX_FOPEN(stream, "symbol.eps", "w");
    solverDraw(pastix_data->solvmatr,
               stream,
               iparm[IPARM_VERBOSE]);
    fclose(stream);
#endif

    /* Let's uncompress the cblk because the solve doesn't know how to deal with compressed information */
    coeftabMemory[spm->flttype-2]( pastix_data->solvmatr );

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_BCSC2CTAB |
                             STEP_SOLVE     |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    iparm[IPARM_START_TASK]++;

    return EXIT_SUCCESS;
}
