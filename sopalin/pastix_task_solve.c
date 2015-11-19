/**
 *
 * @file pastix_task_solve.c
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
#include "order.h"
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

void
pastix_trsm( pastix_data_t *pastix_data,
             int flttype, int side, int uplo, int trans, int diag,
             sopalin_data_t *sopalin_data,
             int nrhs, void *b, int ldb )
{
    switch (flttype) {
    case PastixComplex64:
        sequential_ztrsm( pastix_data, side, uplo, trans, diag,
                             sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sequential_ctrsm( pastix_data, side, uplo, trans, diag,
                             sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sequential_dtrsm( pastix_data, side, uplo, trans, diag,
                             sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sequential_strsm( pastix_data, side, uplo, trans, diag,
                             sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }
}

void
pastix_diag( pastix_data_t *pastix_data, int flttype,
             sopalin_data_t *sopalin_data,
             int nrhs, void *b, int ldb )
{
    switch (flttype) {
    case PastixComplex64:
        sequential_zdiag( pastix_data, sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sequential_cdiag( pastix_data, sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        sequential_ddiag( pastix_data, sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        sequential_sdiag( pastix_data, sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }
}

void
pastix_dsolve( pastix_data_t *pastix_data,
               int flttype, int side, int uplo, int trans, int diag,
               sopalin_data_t *sopalin_data,
               int nrhs, void *b, int ldb )
{
    switch (flttype) {
    case PastixComplex64:
        sequential_z_Dsolve( pastix_data, side, uplo, trans, diag,
                          sopalin_data, nrhs, (pastix_complex64_t *)b, ldb );
        break;
    case PastixComplex32:
        sequential_c_Dsolve( pastix_data, side, uplo, trans, diag,
                          sopalin_data, nrhs, (pastix_complex32_t *)b, ldb );
        break;
    case PastixDouble:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sequential_d_Dsolve( pastix_data, side, uplo, trans, diag,
                          sopalin_data, nrhs, (double *)b, ldb );
        break;
    case PastixFloat:
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sequential_s_Dsolve( pastix_data, side, uplo, trans, diag,
                          sopalin_data, nrhs, (float *)b, ldb );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }
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
pastix_task_solve( pastix_data_t *pastix_data,
                   const pastix_csc_t  *csc,
                   int nrhs, void *b, int ldb )
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

    orderBase( pastix_data->ordemesh, 0 );
    sbackup = solverBackupInit( pastix_data->solvmatr );
    pastix_data->solvmatr->restore = 1;

    switch( pastix_data->bcsc->flttype ) {
    case PastixComplex64:
        z_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->permtab );
        break;

    case PastixComplex32:
        c_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->permtab );
        break;

    case PastixFloat:
        s_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->permtab );
        break;

    case PastixDouble:
    default:
        d_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->permtab );
    }

    {
        sopalin_data_t sopalin_data;
        double timer;

        sopalin_data.solvmtx = pastix_data->solvmatr;

        clockStart(timer);
        switch ( pastix_data->iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLT:
            dump_rhs( "AfterPerm", csc->gN, b );
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans,   PastixNonUnit, &sopalin_data, nrhs, b, ldb );
            dump_rhs( "AfterDown", csc->gN, b );
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixConjTrans, PastixNonUnit, &sopalin_data, nrhs, b, ldb );
            dump_rhs( "AfterUp", csc->gN, b );
            break;

        case PastixFactLDLT:
            dump_rhs( "AfterPerm", csc->gN, b );
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans, PastixUnit, &sopalin_data, nrhs, b, ldb );
            dump_rhs( "AfterDown", csc->gN, b );
            pastix_diag( pastix_data, pastix_data->bcsc->flttype, &sopalin_data, nrhs, b, ldb );
            dump_rhs( "AfterDiag", csc->gN, b );
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixTrans,   PastixUnit, &sopalin_data, nrhs, b, ldb );
            dump_rhs( "AfterUp", csc->gN, b );
            break;

        case PastixFactLDLH:
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans,   PastixUnit, &sopalin_data, nrhs, b, ldb );
            pastix_diag( pastix_data, pastix_data->bcsc->flttype, &sopalin_data, nrhs, b, ldb );
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixConjTrans, PastixUnit, &sopalin_data, nrhs, b, ldb );
            break;

        case PastixFactLU:
        default:
            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans, PastixUnit,    &sopalin_data, nrhs, b, ldb );

            /* pastix_dsolve( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixLower, PastixNoTrans, PastixUnit,    &sopalin_data, nrhs, b, ldb ); */
            /* pastix_dsolve( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit, &sopalin_data, nrhs, b, ldb ); */

            pastix_trsm( pastix_data, pastix_data->bcsc->flttype, PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit, &sopalin_data, nrhs, b, ldb );
            break;
        }
        clockStop(timer);

        pastix_print( 0, 0, OUT_TIME_SOLV, clockVal(timer) );
    }

    switch( pastix_data->bcsc->flttype ) {
    case PastixComplex64:
        z_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->peritab );
        break;

    case PastixComplex32:
        c_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->peritab );
        break;

    case PastixFloat:
        s_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->peritab );
        break;

    case PastixDouble:
    default:
        d_bcscApplyPerm( csc->gN, nrhs, b, ldb,
                         pastix_data->ordemesh->peritab );
    }

    solverBackupRestore( pastix_data->solvmatr, sbackup );
    solverBackupExit( sbackup );

    dump_rhs( "Final", csc->gN, b );

    /* Invalidate following steps, and add factorization step to the ones performed */
    pastix_data->steps &= ~( STEP_SOLVE  |
                             STEP_REFINE );
    pastix_data->steps |= STEP_NUMFACT;

    iparm[IPARM_START_TASK]++;

    return EXIT_SUCCESS;
}
