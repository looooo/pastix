/**
 *
 * @file pastix_task_solve.c
 *
 *  PaStiX solve routines
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
#include "spm.h"
#include "bcsc.h"
#include "order.h"
#include "solver.h"
#include "sopalin_data.h"

#include "z_bcsc.h"
#include "c_bcsc.h"
#include "d_bcsc.h"
#include "s_bcsc.h"

#if defined(PASTIX_DEBUG_SOLVE)
static inline void dump_rhs( char *name, int n, double *b )
{
    int i;
    fprintf(stderr,"%s :", name );
    for (i=0; i<n; i++) {
        if (i%10 == 0)
            fprintf(stderr, "\n");
        fprintf(stderr,"%e ", b[i]);
    }
    fprintf(stderr,"\n");
}
#else
#define dump_rhs(...) do {} while(0)
#endif

int
pastix_subtask_applyorder( pastix_data_t *pastix_data,
                           pastix_coeftype_t flttype, pastix_dir_t dir,
                           pastix_int_t m, pastix_int_t n, void *b, pastix_int_t ldb )
{
    pastix_int_t *perm;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_applyorder: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_applyorder: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        errorPrint("pastix_subtask_applyorder: All steps from pastix_task_init() to pastix_subtask_bcsc2ctab() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Make sure ordering is 0 based */
    orderBase( pastix_data->ordemesh, 0 );

    perm = (dir == PastixDirForward) ? pastix_data->ordemesh->permtab : pastix_data->ordemesh->peritab;

    /* TODO: change name of the ordeing methof since bcsc has nothing to do with this */
    /* See also xlapmr and xlapmt */
    switch( flttype ) {
    case PastixComplex64:
        z_bcscApplyPerm( m, n, b, ldb, perm );
        break;

    case PastixComplex32:
        c_bcscApplyPerm( m, n, b, ldb, perm );
        break;

    case PastixFloat:
        s_bcscApplyPerm( m, n, b, ldb, perm );
        break;

    case PastixDouble:
    default:
        d_bcscApplyPerm( m, n, b, ldb, perm );
    }

    return PASTIX_SUCCESS;
}

int
pastix_subtask_trsm( pastix_data_t *pastix_data,
                     pastix_coeftype_t flttype, pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag,
                     pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    sopalin_data_t sopalin_data;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_trsm: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_trsm: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        errorPrint("pastix_subtask_trsm: All steps from pastix_task_init() to pastix_subtask_bcsc2ctab() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    sopalin_data.solvmtx = pastix_data->solvmatr;

    switch (flttype) {
    case PastixComplex64:
        sopalin_ztrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sopalin_ctrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_dtrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_strsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    return PASTIX_SUCCESS;
}

int
pastix_subtask_diag( pastix_data_t *pastix_data, pastix_coeftype_t flttype,
                     pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
    sopalin_data_t sopalin_data;
    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_subtask_diag: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (b == NULL) {
        errorPrint("pastix_subtask_diag: wrong b parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_BCSC2CTAB) ) {
        errorPrint("pastix_subtask_diag: All steps from pastix_task_init() to pastix_subtask_bcsc2ctab() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    sopalin_data.solvmtx = pastix_data->solvmatr;

    switch (flttype) {
    case PastixComplex64:
        sopalin_zdiag( pastix_data, &sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sopalin_cdiag( pastix_data, &sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        sopalin_ddiag( pastix_data, &sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        sopalin_sdiag( pastix_data, &sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * pastix_task_solve - Solve the given problem.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the solution is stored in place of the right-hand-side vector.
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
pastix_task_solve( pastix_data_t *pastix_data,
                   const pastix_spm_t  *spm,
                   pastix_int_t nrhs, void *b, pastix_int_t ldb )
{
/* #ifdef PASTIX_WITH_MPI */
/*     MPI_Comm       pastix_comm = pastix_data->inter_node_comm; */
/* #endif */
    pastix_int_t  procnum;
    pastix_int_t *iparm;
/*     double        *dparm    = pastix_data->dparm; */
/*     SolverMatrix  *solvmatr = pastix_data->solvmatr; */
    (void)procnum;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_solve: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (spm == NULL) {
        errorPrint("pastix_task_solve: wrong spm parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        errorPrint("pastix_task_solve: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm   = pastix_data->iparm;
    procnum = pastix_data->inter_node_procnum;

    /* Compute P * b */
    pastix_subtask_applyorder( pastix_data, spm->flttype,
                               PastixDirForward, spm->gN, nrhs, b, ldb );

    {
        double timer;

        clockStart(timer);
        switch ( pastix_data->iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLT:
            dump_rhs( "AfterPerm", spm->gN, b );

            /* Solve L y = P b with y = L^t P x */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans,   PastixNonUnit, nrhs, b, ldb );
            dump_rhs( "AfterDown", spm->gN, b );

            /* Solve y = L^t (P x) */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixConjTrans, PastixNonUnit, nrhs, b, ldb );
            dump_rhs( "AfterUp", spm->gN, b );
            break;

        case PastixFactLDLT:
            dump_rhs( "AfterPerm", spm->gN, b );

            /* Solve L y = P b with y = D L^t P x */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans, PastixUnit, nrhs, b, ldb );
            dump_rhs( "AfterDown", spm->gN, b );

            /* Solve y = D z with z = (L^t P x) */
            pastix_subtask_diag( pastix_data, pastix_data->bcsc->flttype, nrhs, b, ldb );
            dump_rhs( "AfterDiag", spm->gN, b );

            /* Solve z = L^t (P x) */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixTrans,   PastixUnit, nrhs, b, ldb );
            dump_rhs( "AfterUp", spm->gN, b );
            break;

        case PastixFactLDLH:
            /* Solve L y = P b with y = D L^h P x */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans,   PastixUnit, nrhs, b, ldb );

            /* Solve y = D z with z = (L^h P x) */
            pastix_subtask_diag( pastix_data, pastix_data->bcsc->flttype, nrhs, b, ldb );

            /* Solve z = L^h (P x) */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixConjTrans, PastixUnit, nrhs, b, ldb );
            break;

        case PastixFactLU:
        default:
            /* Solve L y = P b with y = U P x */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans, PastixUnit,    nrhs, b, ldb );

            /* Solve y = U (P x) */
            pastix_subtask_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit, nrhs, b, ldb );
            break;
        }
        clockStop(timer);

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
            pastix_print( 0, 0, OUT_TIME_SOLV, clockVal(timer) );
        }
    }

    /* Compute P^t * b */
    pastix_subtask_applyorder( pastix_data, spm->flttype,
                               PastixDirBackward, spm->gN, nrhs, b, ldb );
    dump_rhs( "Final", spm->gN, b );

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_SOLVE  |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    return EXIT_SUCCESS;
}
