/**
 *
 * @file pastix_task_solve.c
 *
 *  PaStiX solve routines
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Tony Delarue
 * @date 2022-07-07
 *
 **/
#define _GNU_SOURCE 1
#include "common.h"
#include "bcsc/bcsc.h"
#include "pastix/order.h"
#include "order/order_internal.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"

#include "bcsc/bcsc_z.h"
#include "bcsc/bcsc_c.h"
#include "bcsc/bcsc_d.h"
#include "bcsc/bcsc_s.h"

#include <lapacke.h>

#if defined(PASTIX_DEBUG_SOLVE)
#include <spm/z_spm.h>
#include <spm/c_spm.h>
#include <spm/d_spm.h>
#include <spm/s_spm.h>

static volatile int32_t step = 0;

static inline void
dump_rhs( pastix_data_t  *pastix_data,
          const char     *name,
          spm_coeftype_t  flttype,
          int             m,
          int             n,
          const void     *b,
          int             ldb )
{
    FILE *f;
    f = pastix_fopenw( (pastix_data->dir_local == NULL) ? "." : pastix_data->dir_local,
                       name, "w" );

    switch( flttype ) {
    case SpmComplex64:
        z_spmDensePrint( f, m, n, b, ldb );
        break;
    case SpmComplex32:
        c_spmDensePrint( f, m, n, b, ldb );
        break;
    case SpmDouble:
        d_spmDensePrint( f, m, n, b, ldb );
        break;
    case SpmFloat:
        s_spmDensePrint( f, m, n, b, ldb );
        break;
    case SpmPattern:
        ;
    }

    fclose(f);
}
#else
static inline void
dump_rhs( pastix_data_t  *pastix_data,
          const char     *name,
          spm_coeftype_t  flttype,
          int             m,
          int             n,
          const void     *b,
          int             ldb )
{
    (void)pastix_data;
    (void)name;
    (void)m;
    (void)n;
    (void)b;
    (void)flttype;
    (void)ldb;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Initialize an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] B_ptr
 *          On entry, an allocated pastix_rhs_t data structure.
 *          On exit, the data is initialized to be used by the pastix_subtask_*
 *          functions.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsInit( pastix_data_t *pastix_data,
               pastix_rhs_t  *B_ptr )
{
    pastix_rhs_t B;

    if ( pastix_data == NULL ) {
        pastix_print_error( "pastixRhsInit: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( B_ptr == NULL ) {
        pastix_print_error( "pastixRhsInit: wrong B parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    *B_ptr = malloc( sizeof(struct pastix_rhs_s) );
    B = *B_ptr;

    B->allocated = 0;
    B->flttype   = PastixPattern;
    B->m         = -1;
    B->n         = -1;
    B->ld        = -1;
    B->b         = NULL;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Cleanup an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] B
 *          On entry, the initialized pastix_rhs_t data structure.
 *          On exit, the structure is destroyed and should no longer be used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsFinalize( pastix_data_t *pastix_data,
                   pastix_rhs_t   B )
{
    if ( pastix_data == NULL ) {
        pastix_print_error( "pastixRhsFinalize: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( B == NULL ) {
        pastix_print_error( "pastixRhsFinalize: wrong B parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( B->b != NULL ) {
        if ( B->allocated ) {
            free( B->b );
        }
        else {
            pastix_print_warning( "Calling pastixRhsFinalize before restoring the ordering of vector b.\n"
                                  "Please call:\n"
                                  "  pastix_subtask_applyorder( pastix_data, flttype, PastixDirBackward, m, n,\n"
                                  "                             b, ldb, Bp );\n"
                                  "prior to this call to restore it.\n" );
        }
    }

    free( B );
    return PASTIX_SUCCESS;
}

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
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
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
pastix_subtask_applyorder( pastix_data_t    *pastix_data,
                           pastix_coeftype_t flttype,
                           pastix_dir_t      dir,
                           pastix_int_t      m,
                           pastix_int_t      n,
                           void             *b,
                           pastix_int_t      ldb,
                           pastix_rhs_t      Bp )
{
    pastix_int_t *perm = NULL;
    int ts;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_subtask_applyorder: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (Bp == NULL) {
        pastix_print_error( "pastix_subtask_applyorder: wrong Bp parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_CSC2BCSC) ) {
        pastix_print_error( "pastix_subtask_applyorder: All steps from pastix_task_init() to pastix_subtask_csc2bcsc() have to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Make sure ordering is 0 based */
    if ( pastix_data->ordemesh->baseval != 0 ) {
        pastix_print_error( "pastix_subtask_applyorder: ordermesh must be 0-based" );
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( dir == PastixDirForward ) {
        Bp->flttype = flttype;
        Bp->m       = m;
        Bp->n       = n;
        Bp->ld      = Bp->m;
        Bp->b       = b;
    }
    else {
        assert( Bp->flttype == flttype );
        assert( Bp->b  == b );
        assert( Bp->m  == m );
        assert( Bp->n  == n );
        assert( Bp->ld >= Bp->m );
    }

    ts   = pastix_data->iparm[IPARM_APPLYPERM_WS];
    perm = orderGetExpandedPeritab( pastix_data->ordemesh, pastix_data->csc );

    /* See also xlapmr and xlapmt */
    switch( flttype ) {
    case PastixComplex64:
        bvec_zlapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixComplex32:
        bvec_clapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixFloat:
        bvec_slapmr( ts, dir, m, n, b, ldb, perm );
        break;

    case PastixDouble:
    default:
        bvec_dlapmr( ts, dir, m, n, b, ldb, perm );
    }

#if defined(PASTIX_DEBUG_SOLVE)
    {
        char *fullname;
        int32_t lstep = pastix_atomic_inc_32b( &step );
        asprintf( &fullname, "solve_%02d_applyorder_%s.rhs", lstep,
                  ( dir == PastixDirForward ) ? "Forward" : "Backward" );

        dump_rhs( pastix_data, fullname,
                  flttype, m, n, b, ldb );
        free(fullname);
    }
#endif

    if ( dir == PastixDirBackward ) {
        if ( Bp->allocated ) {
            free( Bp->b );
        }

        Bp->allocated = 0;
        Bp->flttype   = PastixPattern;
        Bp->m         = -1;
        Bp->n         = -1;
        Bp->ld        = -1;
        Bp->b         = NULL;
    }

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
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
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
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vector (can be multiple RHS).
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
pastix_subtask_trsm( pastix_data_t    *pastix_data,
                     pastix_side_t     side,
                     pastix_uplo_t     uplo,
                     pastix_trans_t    trans,
                     pastix_diag_t     diag,
                     pastix_rhs_t      Bp )
{
    sopalin_data_t    sopalin_data;
    pastix_int_t      nrhs, i, bs, ldb;
    pastix_coeftype_t flttype;
    void             *b;

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
    nrhs = Bp->n;
    b    = Bp->b;
    ldb  = Bp->ld;
    bs   = nrhs;
#if defined(PASTIX_WITH_MPI)
    bs = 1;
#endif

    /*
     * Ensure that the scheduler is correct and is in the same
     * family that the one used for the factorization.
     */
    pastix_check_and_correct_scheduler( pastix_data );

    sopalin_data.solvmtx = pastix_data->solvmatr;

    switch (flttype) {
    case PastixComplex64:
    {
        pastix_complex64_t *lb = b;
        for( i = 0; i < nrhs; i+=bs, lb += ldb ) {
            sopalin_ztrsm( pastix_data, side, uplo, trans, diag,
                           &sopalin_data, bs, lb, ldb );
        }
    }
    break;
    case PastixComplex32:
    {
        pastix_complex32_t *lb = b;
        for( i = 0; i < nrhs; i+=bs, lb += ldb ) {
            sopalin_ctrsm( pastix_data, side, uplo, trans, diag,
                           &sopalin_data, bs, lb, ldb );
        }
    }
    break;
    case PastixDouble:
    {
        double *lb = b;
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        for( i = 0; i < nrhs; i+=bs, lb += ldb ) {
            sopalin_dtrsm( pastix_data, side, uplo, trans, diag,
                           &sopalin_data, bs, lb, ldb );
        }
    }
    break;
    case PastixFloat:
    {
        float *lb = b;
        trans = (trans == PastixConjTrans) ? PastixTrans : trans;
        for( i = 0; i < nrhs; i+=bs, lb += ldb ) {
            sopalin_strsm( pastix_data, side, uplo, trans, diag,
                           &sopalin_data, bs, lb, ldb );
        }
    }
    break;
    default:
        fprintf(stderr, "Unknown floating point arithmetic\n" );
    }

#if defined(PASTIX_DEBUG_SOLVE)
    {
        char *fullname;
        int32_t lstep = pastix_atomic_inc_32b( &step );
        asprintf( &fullname, "solve_%02d_trsm_%c%c%c%c.rhs", lstep,
                  (side  == PastixLeft)  ? 'L' : 'R',
                  (uplo  == PastixUpper) ? 'U' : 'L',
                  (trans == PastixConjTrans) ? 'C' : (trans == PastixTrans ? 'T' : 'N'),
                  (diag  == PastixUnit)  ? 'U' : 'N' );

        dump_rhs( pastix_data, fullname,
                  flttype, pastix_data->bcsc->gN, nrhs, b, ldb );
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
 * @param[in] flttype
 *          This arithmetic of the sparse matrix.
 *
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vector (can be multiple RHS).
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
pastix_subtask_diag( pastix_data_t *pastix_data,
                     pastix_rhs_t   Bp )
{
    sopalin_data_t sopalin_data;
    pastix_coeftype_t flttype;
    pastix_int_t      nrhs;
    void             *b;
    pastix_int_t      ldb;

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

    flttype = Bp->flttype;
    nrhs = Bp->n;
    b    = Bp->b;
    ldb  = Bp->ld;

    /*
     * Ensure that the scheduler is correct and is in the same
     * family that the one used for the factorization.
     */
    pastix_check_and_correct_scheduler( pastix_data );

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

#if defined(PASTIX_DEBUG_SOLVE)
    {
        char *fullname;
        int32_t lstep = pastix_atomic_inc_32b( &step );
        asprintf( &fullname, "solve_%02d_diag.rhs", lstep );

        dump_rhs( pastix_data, fullname,
                  flttype, pastix_data->bcsc->gN, nrhs, b, ldb );
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
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
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
        double timer;
        /* Start timer */
        clockSyncStart( timer, pastix_data->inter_node_comm );

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

        pastix_subtask_trsm( pastix_data, PastixLeft, uplo, trans, diag, Bp );

        /*
         * Solve the diagonal step
         */
        if( (factotype == PastixFactLDLT) ||
            (factotype == PastixFactLDLH) )
        {
            /* Solve y = D z with z = ([L^t | L^h] P x) */
            pastix_subtask_diag( pastix_data, Bp );
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

        pastix_subtask_trsm( pastix_data,
                             PastixLeft, uplo, trans, diag, Bp );

        /* Stop Timer */
        clockSyncStop( timer, pastix_data->inter_node_comm );

        pastix_data->dparm[DPARM_SOLV_TIME] = clockVal(timer);
        if ( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            pastix_print( pastix_data->inter_node_procnum, 0, OUT_TIME_SOLV,
                          pastix_data->dparm[DPARM_SOLV_TIME] );
        }
    }

    return EXIT_SUCCESS;
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
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
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
 * @param[in] nrhs
 *          The number of right-and-side vectors.
 *
 * @param[inout] b
 *          The right-and-side vectors (can be multiple RHS).
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
                   pastix_int_t   nrhs,
                   void          *b,
                   pastix_int_t   ldb )
{
    pastix_bcsc_t *bcsc;
    void          *bglob = NULL;
    void          *tmp   = NULL;
    void          *sb;
    pastix_rhs_t   Bp;

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

    bcsc = pastix_data->bcsc;

    /* The spm is distributed, so we have to gather the RHS */
    if ( pastix_data->csc->loc2glob != NULL ) {
        if( pastix_data->iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
            pastix_print( pastix_data->procnum, 0, "pastix_task_solve: the RHS has to be centralized for the moment\n" );
        }

        tmp = b;
        bglob = malloc( pastix_data->csc->gNexp * nrhs * pastix_size_of( bcsc->flttype ) );
        spmGatherRHS( nrhs, pastix_data->csc, b, ldb,
                      -1, bglob, pastix_data->csc->gNexp );
        b   = bglob;
        ldb = pastix_data->csc->gNexp;
    }

    /* Generating halved-precision vector */
    if ( pastix_data->iparm[IPARM_MIXED] &&
         ((bcsc->flttype == PastixComplex64) || (bcsc->flttype == PastixDouble)) )
    {
        int    rc;
        size_t size = ldb * nrhs;

        switch(bcsc->flttype) {
        case PastixComplex64:
            sb = malloc( size * sizeof(pastix_complex32_t) );
            rc = LAPACKE_zlag2c_work( LAPACK_COL_MAJOR, ldb, nrhs,
                                      b, ldb, sb, ldb );
            break;
        case PastixDouble:
            sb = malloc( size * sizeof(float) );
            rc = LAPACKE_dlag2s_work( LAPACK_COL_MAJOR, ldb, nrhs,
                                      b, ldb, sb, ldb );
            break;
        default:
            rc = 1;
            pastix_print_error( "pastix_task_solve: Invalid float type for mixed-precision" );
        }

        if ( rc ) {
            free( sb );
            return PASTIX_ERR_INTERNAL;
        }
    }
    else {
        sb = b;
    }

    /* Compute P * b */
    pastixRhsInit( pastix_data, &Bp );
    pastix_subtask_applyorder( pastix_data, pastix_data->solvmatr->flttype,
                               PastixDirForward, bcsc->gN, nrhs, sb, ldb, Bp );

    /* Solve A x = b */
    pastix_subtask_solve( pastix_data, Bp );

    /* Compute P^t * b */
    pastix_subtask_applyorder( pastix_data, pastix_data->solvmatr->flttype,
                               PastixDirBackward, bcsc->gN, nrhs, sb, ldb, Bp );
    pastixRhsFinalize( pastix_data, Bp );

    /* Freeing mixed-precision vectors and reverting to given precision */
    if ( pastix_data->iparm[IPARM_MIXED] &&
         ((bcsc->flttype == PastixComplex64) || (bcsc->flttype == PastixDouble)) )
    {
        int rc;

        switch(bcsc->flttype) {
        case PastixComplex64:
            rc = LAPACKE_clag2z_work( LAPACK_COL_MAJOR, ldb, nrhs,
                                      sb, ldb, b, ldb );
            break;
        case PastixDouble:
            rc = LAPACKE_slag2d_work( LAPACK_COL_MAJOR, ldb, nrhs,
                                      sb, ldb, b, ldb );
            break;
        default:
            pastix_print_error( "pastix_task_solve: Invalid float type for mixed-precision" );
        }

        assert( b != sb );
        free( sb );
        if( rc ) {
            return PASTIX_ERR_INTERNAL;
        }
    }

    if( tmp != NULL ) {
        pastix_int_t ldbglob = ldb;
        ldb = pastix_data->csc->nexp;
        b   = tmp;

        spmExtractLocalRHS( nrhs, pastix_data->csc, bglob, ldbglob, b, ldb );
        memFree_null( bglob );
    }

    return EXIT_SUCCESS;
}
