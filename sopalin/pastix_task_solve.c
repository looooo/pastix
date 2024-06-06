/**
 *
 * @file pastix_task_solve.c
 *
 *  PaStiX solve routines
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @author Tony Delarue
 * @date 2023-11-10
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#include "bcsc/bcsc.h"
#include "pastix/order.h"
#include "order/order_internal.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_papi.h"

#include "bcsc/bcsc_z.h"
#include "bcsc/bcsc_c.h"
#include "bcsc/bcsc_d.h"
#include "bcsc/bcsc_s.h"

#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#if defined(PASTIX_DEBUG_SOLVE)
#include <spm/z_spm.h>
#include <spm/c_spm.h>
#include <spm/d_spm.h>
#include <spm/s_spm.h>

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] name
 *          TODO
 *
 * @param[in] B
 *          TODO
 *
 *******************************************************************************/
static inline void
pastix_rhs_dump( pastix_data_t *pastix_data,
                 const char    *name,
                 pastix_rhs_t   B )
{
    static volatile int32_t step = 0;

    int32_t lstep = pastix_atomic_inc_32b( &step );
    char   *fullname;
    FILE   *f;
    int     rc;

    rc = asprintf( &fullname, "%02d_%s.%02d.rhs",
                   lstep, name, pastix_data->solvmatr->clustnum );
    assert( rc != -1 );

    f = pastix_fopenw( (pastix_data->dir_local == NULL) ? "." : pastix_data->dir_local,
                       fullname, "w" );
    free( fullname );

    switch( B->flttype ) {
    case SpmComplex64:
        z_spmDensePrint( f, B->m, B->n, B->b, B->ld );
        break;
    case SpmComplex32:
        c_spmDensePrint( f, B->m, B->n, B->b, B->ld );
        break;
    case SpmDouble:
        d_spmDensePrint( f, B->m, B->n, B->b, B->ld );
        break;
    case SpmFloat:
        s_spmDensePrint( f, B->m, B->n, B->b, B->ld );
        break;
    case SpmPattern:
        ;
    }

    fclose(f);
    (void)rc;
}
#else
static inline void
pastix_rhs_dump( pastix_data_t *pastix_data,
                 const char    *name,
                 pastix_rhs_t   B )
{
    (void)pastix_data;
    (void)name;
    (void)B;
}
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a permutation on the right-and-side vector before the solve step.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_APPLYPERM_WS.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] dir
 *          Forward or backword application of the permutation.
 *
 * @param[in] m
 *          Size of the right-and-side vectors.
 *
 * @param[in] n
 *          Number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors of size ldb-by-n.
 *          If dir == PastixDirForward, b is used as input to initialize the
 *          permuted b. If the matrix is replicated on all nodes or in shared
 *          memory, the permuted vector has the same size as b, and b is
 *          modified in-place.
 *          If dir == PastixDirBackward, b must be allocated on entry and is
 *          filled by the reverse permutation of Bp. Note that if the matrix is
 *          replicated on all nodes or in shared memory, b is modified in-place.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 * @param[inout] Bp
 *          The right-and-side vectors of size ldb-by-n.  On entry, must be
 *          initialized through pastixRhsInit(). On exit, contains the
 *          information about the permuted vector when dir ==
 *          PastixDirForward. When dir == PastixDirBackward, the data is used as
 *          input to update b.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_applyorder( pastix_data_t *pastix_data,
                           pastix_dir_t   dir,
                           pastix_int_t   m,
                           pastix_int_t   n,
                           void          *b,
                           pastix_int_t   ldb,
                           pastix_rhs_t   Bp )
{
    pastix_coeftype_t flttype;

    /*
     * Checks parameters.
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_subtask_applyorder: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (Bp == NULL) {
        pastix_print_error( "pastix_subtask_applyorder: wrong Bp parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Makes sure ordering is 0 based. */
    if ( pastix_data->ordemesh->baseval != 0 ) {
        pastix_print_error( "pastix_subtask_applyorder: ordermesh must be 0-based" );
        return PASTIX_ERR_BADPARAMETER;
    }

    flttype = pastix_data->csc->flttype;
#if defined(PASTIX_DEBUG_SOLVE)
    if ( dir == PastixDirForward ) {
        struct pastix_rhs_s rhsb = {
            .allocated = 0,
            .flttype   = flttype,
            .m         = m,
            .n         = n,
            .b         = b,
            .ld        = ldb,
        };
        pastix_rhs_dump( pastix_data, "applyorder_Forward_input", &rhsb );
    }
    else {
        pastix_rhs_dump( pastix_data, "applyorder_Backward_input", Bp );
    }
#endif

    /* See also xlapmr and xlapmt */
    switch( flttype ) {
    case PastixComplex64:
        bvec_zlapmr( pastix_data, dir, m, n, b, ldb, Bp );
        break;

    case PastixComplex32:
        bvec_clapmr( pastix_data, dir, m, n, b, ldb, Bp );
        break;

    case PastixDouble:
        bvec_dlapmr( pastix_data, dir, m, n, b, ldb, Bp );
        break;

    case PastixFloat:
        bvec_slapmr( pastix_data, dir, m, n, b, ldb, Bp );
        break;

    default:
        pastix_print_error( "The floating type of the rhs is not defined\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

#if defined(PASTIX_DEBUG_SOLVE)
    if ( dir == PastixDirForward ) {
        pastix_rhs_dump( pastix_data, "applyorder_Forward_output", Bp );
    }
    else {
        struct pastix_rhs_s rhsb = {
            .allocated = 0,
            .flttype   = flttype,
            .m         = m,
            .n         = n,
            .b         = b,
            .ld        = ldb,
        };
        pastix_rhs_dump( pastix_data, "applyorder_Backward_output", &rhsb );
    }
#endif

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a triangular solve on the right-and-side vectors.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] side
 *          Left or right application.
 *
 * @param[in] uplo
 *          Upper or Lower part.
 *
 * @param[in] trans
 *          With or without transposition (or conjugate transposition).
 *
 * @param[in] diag
 *          Diagonal terms are unit or not.
 *
 * @param[inout] Bp
 *          The right-and-side vector (can be multiple rhs).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_trsm( pastix_data_t *pastix_data,
                     pastix_side_t  side,
                     pastix_uplo_t  uplo,
                     pastix_trans_t trans,
                     pastix_diag_t  diag,
                     pastix_rhs_t   Bp )
{
    SolverMatrix     *solvmtx;
    sopalin_data_t    sopalin_data;
    pastix_coeftype_t flttype;
    int               prevnt;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_subtask_trsm: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (Bp == NULL) {
        pastix_print_error( "pastix_subtask_trsm: wrong Bp parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        pastix_print_error( "pastix_subtask_trsm: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    flttype = Bp->flttype;
    solvmtx = pastix_data->solvmatr;

    if ( Bp->cblkb == NULL ) {
        pastix_int_t nbbuffers;

        nbbuffers = solvmtx->faninnbr + solvmtx->recvnbr;
        if ( nbbuffers > 0 ) {
            Bp->cblkb = calloc( nbbuffers, sizeof(void*) );
        }
    }

    /*
     * Ensure that the scheduler is correct and is in the same
     * family that the one used for the factorization.
     */
    pastix_check_and_correct_scheduler( pastix_data );

    sopalin_data.solvmtx = solvmtx;

    prevnt = pastixBlasSetNumThreadsOne();

    switch (flttype) {
    case PastixComplex64:
    {
        sopalin_ztrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, Bp );
    }
    break;
    case PastixComplex32:
    {
        sopalin_ctrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, Bp );
    }
    break;
    case PastixDouble:
    {
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_dtrsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, Bp );
    }
    break;
    case PastixFloat:
    {
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        sopalin_strsm( pastix_data, side, uplo, trans, diag,
                       &sopalin_data, Bp );
    }
    break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    pastixBlasSetNumThreads( prevnt );

#if defined(PASTIX_DEBUG_SOLVE)
    {
        char *fullname;
        asprintf( &fullname, "solve_trsm_%c%c%c%c",
                  (side  == PastixLeft)  ? 'L' : 'R',
                  (uplo  == PastixUpper) ? 'U' : 'L',
                  (trans == PastixConjTrans) ? 'C' : (trans == PastixTrans ? 'T' : 'N'),
                  (diag  == PastixUnit)  ? 'U' : 'N' );

        pastix_rhs_dump( pastix_data, fullname, Bp );
        free(fullname);
    }
#endif

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Apply a diagonal operation on the right-and-side vectors.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] Bp
 *          The right-and-side vector (can be multiple rhs).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_diag( pastix_data_t *pastix_data,
                     pastix_rhs_t   Bp )
{
    sopalin_data_t sopalin_data;
    int            prevnt;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_subtask_diag: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (Bp == NULL) {
        pastix_print_error( "pastix_subtask_diag: wrong Bp parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        pastix_print_error( "pastix_subtask_trsm: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Ensure that the scheduler is correct and is in the same
     * family that the one used for the factorization.
     */
    pastix_check_and_correct_scheduler( pastix_data );

    sopalin_data.solvmtx = pastix_data->solvmatr;

    prevnt = pastixBlasSetNumThreadsOne();

    switch ( Bp->flttype ) {
    case PastixComplex64:
        sopalin_zdiag( pastix_data, &sopalin_data, Bp );
        break;
    case PastixComplex32:
        sopalin_cdiag( pastix_data, &sopalin_data, Bp );
        break;
    case PastixDouble:
        sopalin_ddiag( pastix_data, &sopalin_data, Bp );
        break;
    case PastixFloat:
        sopalin_sdiag( pastix_data, &sopalin_data, Bp );
        break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

    pastixBlasSetNumThreads( prevnt );

    pastix_rhs_dump( pastix_data, "solve_diag", Bp );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Solve the given problem without applying the permutation.
 *
 * @warning The input vector is considered already permuted. For a solve step
 * with permutation, see pastix_task_solve()
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] transA
 *          PastixNoTrans:   A is not transposed (CSC matrix)
 *          PastixTrans:     A is transposed (CSR of symmetric/general matrix)
 *          PastixConjTrans: A is conjugate transposed (CSR of hermitian matrix)
 *
 * @param[inout] Bp
 *          The right-and-side vectors (can be multiple rhs).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_solve_adv( pastix_data_t  *pastix_data,
                          pastix_trans_t  transA,
                          pastix_rhs_t    Bp )
{
    pastix_bcsc_t     *bcsc;
    pastix_factotype_t factotype;

    pastix_uplo_t  uplo;
    pastix_diag_t  diag;
    pastix_trans_t trans;
    pastix_trans_t transfact = PastixTrans;
    pastix_rhs_t   sBp;
    pastix_rhs_t   B;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_task_solve: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        pastix_print_error( "pastix_task_solve: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    bcsc      = pastix_data->bcsc;
    factotype = pastix_data->iparm[IPARM_FACTORIZATION];

    /**
     *
     * Summary of the operations to perform the full solve if the original A was
     * either transposed or conjugate transposed.
     *
     * op(A)     | Factorization       | Step 1    | Step 2
     * ----------+---------------------+-----------+-----------
     * NoTrans   | L U                 | L   y = b | U   x = y
     * NoTrans   | L L^t               | L   y = b | L^t x = y
     * NoTrans   | L L^h               | L   y = b | L^h x = y
     * Trans     |(L U  )^t = U^t  L^t | U^t y = b | L^t x = y
     * Trans     |(L L^t)^t = L    L^t | L   y = b | L^t x = y
     * Trans     |(L L^h)^t = c(L) L^t | Not handled (c(L))
     * ConjTrans |(L U  )^h = U^h  L^h | Not handled (U^h)
     * ConjTrans |(L L^t)^h = c(L) L^h | Not handled (c(L))
     * ConjTrans |(L L^h)^h = L    L^h | L   y = b | L^h x = y
     *
     */
    /* Value of the transpose case */
    if ( ((bcsc->flttype == PastixComplex32) || (bcsc->flttype == PastixComplex64)) &&
         ((factotype == PastixFactLLH ) || ( factotype == PastixFactLDLH )) )
    {
        transfact = PastixConjTrans;
    }

    if ( (transA != PastixNoTrans) &&
         (transA != transfact) )
    {
        pastix_print_error( "pastix_task_solve: transA incompatible with the factotype used (require extra conj(L) not handled)" );
        return PASTIX_ERR_BADPARAMETER;
    }

    {
        double timer, energy;

        papiEnergyStart();
        if ( pastix_data->iparm[IPARM_TRACE] & PastixTraceSolve ) {
            kernelsTraceStart();
        }
        /* Start timer */
        clockSyncStart( timer, pastix_data->inter_node_comm );

        if ( pastix_data->iparm[IPARM_MIXED] &&
             ( ( Bp->flttype == PastixComplex64 ) || ( Bp->flttype == PastixDouble ) ) )
        {
            pastixRhsInit( &sBp );
            pastixRhsDoubletoSingle( Bp, sBp );

            B = sBp;
        }
        else {
            B = Bp;
        }

        /*
         * Solve the first step
         */
        uplo  = PastixLower;
        trans = PastixNoTrans;
        diag  = PastixNonUnit;

        if (( transA != PastixNoTrans ) && ( factotype == PastixFactLU )) {
            uplo  = PastixUpper;
            trans = transA;
        }
        if (( transA == PastixNoTrans ) && ( factotype == PastixFactLU )) {
            diag = PastixUnit;
        }

        if (( factotype == PastixFactLDLT ) ||
            ( factotype == PastixFactLDLH ) )
        {
            diag = PastixUnit;
        }

        pastix_subtask_trsm( pastix_data, PastixLeft, uplo, trans, diag, B );

        /*
         * Solve the diagonal step
         */
        if( (factotype == PastixFactLDLT) ||
            (factotype == PastixFactLDLH) )
        {
            /* Solve y = D z with z = ([L^t | L^h] P x) */
            pastix_subtask_diag( pastix_data, B );
        }

        /*
         * Solve the second step
         */
        uplo  = PastixLower;
        trans = transfact;
        diag  = PastixNonUnit;

        if (( transA == PastixNoTrans ) && ( factotype == PastixFactLU )) {
            uplo  = PastixUpper;
            trans = PastixNoTrans;
        }
        if (( transA != PastixNoTrans ) && ( factotype == PastixFactLU )) {
            diag = PastixUnit;
        }

        if (( factotype == PastixFactLDLT ) ||
            ( factotype == PastixFactLDLH ) )
        {
            diag = PastixUnit;
        }

        pastix_subtask_trsm( pastix_data, PastixLeft, uplo, trans, diag, B );

        if ( pastix_data->iparm[IPARM_MIXED] &&
             ( ( Bp->flttype == PastixComplex64 ) || ( Bp->flttype == PastixDouble ) ) )
        {
            pastixRhsSingleToDouble( sBp, Bp );
            pastixRhsFinalize( sBp );
        }

        /* Stop Timer */
        clockSyncStop( timer, pastix_data->inter_node_comm );
        if ( pastix_data->iparm[IPARM_TRACE] & PastixTraceSolve ) {
            kernelsTraceStop();
        }
        energy = papiEnergyStop();

        pastix_data->dparm[DPARM_SOLV_TIME]   = clockVal(timer);
        pastix_data->dparm[DPARM_SOLV_ENERGY] = energy;

        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            pastix_print( pastix_data->inter_node_procnum, 0, OUT_TIME_SOLV,
                          pastix_data->dparm[DPARM_SOLV_TIME] );
#if defined(PASTIX_WITH_PAPI)
            pastix_print( pastix_data->inter_node_procnum, 0, OUT_SOLVE_ENERGY,
                          pastix_print_value_deci( pastix_data->dparm[DPARM_SOLV_ENERGY] ),
                          pastix_print_unit_deci(  pastix_data->dparm[DPARM_SOLV_ENERGY] ),
                          pastix_print_value_deci( pastix_data->dparm[DPARM_SOLV_ENERGY] / pastix_data->dparm[DPARM_SOLV_TIME] ),
                          pastix_print_unit_deci(  pastix_data->dparm[DPARM_SOLV_ENERGY] / pastix_data->dparm[DPARM_SOLV_TIME] ) );
#endif
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Solve the given problem without applying the permutation.
 *
 * @warning The input vector is considered already permuted. For a solve step
 * with permutation, see pastix_task_solve()
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_TRANSPOSE_SOLVE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] Bp
 *          The right-and-side vectors (can be multiple rhs).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_subtask_solve( pastix_data_t *pastix_data,
                      pastix_rhs_t   Bp )
{
    return pastix_subtask_solve_adv( pastix_data,
                                     pastix_data->iparm[IPARM_TRANSPOSE_SOLVE],
                                     Bp );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Solve the given problem.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_FACTORIZATION, IPARM_TRANSPOSE_SOLVE.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[in] m
 *          The local number of rows of the right-and-side vectors.
 *          If the spm is centralized or replicated, m must be equal to
 *          spm->gNexp, otherwise m must be equal to spm->nexp.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple rhs).
 *          On exit, the solution is stored in place of the right-hand-side vector.
 *
 * @param[in] ldb
 *          The leading dimension of the right-and-side vectors.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastix_task_solve( pastix_data_t *pastix_data,
                   pastix_int_t   m,
                   pastix_int_t   nrhs,
                   void          *b,
                   pastix_int_t   ldb )
{
    pastix_rhs_t Bp;
    pastix_int_t rc;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_task_solve: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        pastix_print_error( "pastix_task_solve: Numerical factorization hasn't been done." );
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Compute P * b */
    rc = pastixRhsInit( &Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    rc = pastix_subtask_applyorder( pastix_data, PastixDirForward, m, nrhs, b, ldb, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Solve A x = b */
    rc = pastix_subtask_solve( pastix_data, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Compute P^t * b */
    rc = pastix_subtask_applyorder( pastix_data, PastixDirBackward, m, nrhs, b, ldb, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    rc = pastixRhsFinalize( Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    return rc;
}
