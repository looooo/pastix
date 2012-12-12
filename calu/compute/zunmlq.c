/**
 *
 * @file zunmlq.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zunmlq - Overwrites the general complex M-by-N matrix C with
 * 
 *                  SIDE = 'L'     SIDE = 'R'
 *  TRANS = 'N':      Q * C          C * Q
 *  TRANS = 'C':      Q**H * C       C * Q**H
 * 
 *  where Q is a complex unitary matrix defined as the product of k
 *  elementary reflectors
 * 
 *        Q = H(1) H(2) . . . H(k)
 * 
 *  as returned by PLASMA_zgeqrf. Q is of order M if SIDE = PlasmaLeft
 *  and of order N if SIDE = PlasmaRight.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = PlasmaLeft:  apply Q or Q**H from the left;
 *          = PlasmaRight: apply Q or Q**H from the right.
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   no transpose, apply Q;
 *          = PlasmaConjTrans: conjugate transpose, apply Q**H.
 *
 * @param[in] M
 *          The number of rows of the matrix C. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix C. N >= 0.
 *
 * @param[in] K
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          If side == PlasmaLeft,  M >= K >= 0.
 *          If side == PlasmaRight, N >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K).
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by PLASMA_zgelqf.
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q**H*C.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zunmlq_Tile
 * @sa PLASMA_zunmlq_Tile_Async
 * @sa PLASMA_cunmlq
 * @sa PLASMA_dormlq
 * @sa PLASMA_sormlq
 * @sa PLASMA_zgelqf
 *
 ******************************************************************************/
int PLASMA_zunmlq(PLASMA_enum side, PLASMA_enum trans, int M, int N, int K,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_desc *descT,
                  PLASMA_Complex64_t *C, int LDC)
{
    int NB, An;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descC;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zunmlq", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if (side == PlasmaLeft)
        An = M;
    else 
        An = N;

    /* Check input arguments */
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        plasma_error("PLASMA_zunmlq", "illegal value of side");
        return -1;
    }
    if ((trans != PlasmaConjTrans) && (trans != PlasmaNoTrans)){
        plasma_error("PLASMA_zunmlq", "illegal value of trans");
        return -2;
    }
    if (M < 0) {
        plasma_error("PLASMA_zunmlq", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_zunmlq", "illegal value of N");
        return -4;
    }
    if ((K < 0) || (K > An)) {
        plasma_error("PLASMA_zunmlq", "illegal value of K");
        return -5;
    }
    if (LDA < max(1, K)) {
        plasma_error("PLASMA_zunmlq", "illegal value of LDA");
        return -7;
    }
    if (LDC < max(1, M)) {
        plasma_error("PLASMA_zunmlq", "illegal value of LDC");
        return -10;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGELS, M, K, N);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zunmlq", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT & NTRHS */
    NB   = PLASMA_NB;
    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, K, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descC, C, NB, NB, LDC, N, 0, 0, M,  N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, An, 0, 0, K, An,
                            sequence, &request);
        plasma_ziplap2tile( descC, C, NB, NB, LDC, N, 0, 0, M,  N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zunmlq_Tile_Async(
        side, trans, &descA, descT, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, An, sequence, &request);
        plasma_ziptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zunmlq_Tile - overwrites the general M-by-N matrix C with Q*C, where Q is an orthogonal
 *  matrix (unitary in the complex case) defined as the product of elementary reflectors returned
 *  by PLASMA_zgelqf_Tile Q is of order M.
 *  All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Intended usage:
 *          = PlasmaLeft:  apply Q or Q**H from the left;
 *          = PlasmaRight: apply Q or Q**H from the right.
 *          Currently only PlasmaLeft is supported.
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   no transpose, apply Q;
 *          = PlasmaConjTrans: conjugate transpose, apply Q**H.
 *          Currently only PlasmaConjTrans is supported.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by PLASMA_zgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_zgelqf.
 *
 * @param[in,out] C
 *          On entry, the M-by-N matrix C.
 *          On exit, C is overwritten by Q*C or Q**H*C.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zunmlq
 * @sa PLASMA_zunmlq_Tile_Async
 * @sa PLASMA_cunmlq_Tile
 * @sa PLASMA_dormlq_Tile
 * @sa PLASMA_sormlq_Tile
 * @sa PLASMA_zgelqf_Tile
 *
 ******************************************************************************/
int PLASMA_zunmlq_Tile(PLASMA_enum side, PLASMA_enum trans,
                       PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zunmlq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zunmlq_Tile_Async(side, trans, A, T, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  Non-blocking equivalent of PLASMA_zunmlq_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_zunmlq
 * @sa PLASMA_zunmlq_Tile
 * @sa PLASMA_cunmlq_Tile_Async
 * @sa PLASMA_dormlq_Tile_Async
 * @sa PLASMA_sormlq_Tile_Async
 * @sa PLASMA_zgelqf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zunmlq_Tile_Async(PLASMA_enum side, PLASMA_enum trans,
                             PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *C,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descC;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zunmlq_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zunmlq_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zunmlq_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zunmlq_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zunmlq_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zunmlq_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descC = *C;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descC.nb != descC.mb) {
        plasma_error("PLASMA_zunmlq_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((side != PlasmaLeft) && (side != PlasmaRight)) {
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((trans != PlasmaConjTrans) && (trans != PlasmaNoTrans)){
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, C, LDC ) */
/*
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        if ( (trans == PlasmaConjTrans) &&
             (side == PlasmaLeft) ) {
            plasma_parallel_call_7(plasma_pzunmlq,
                PLASMA_enum, side,
                PLASMA_enum, trans,
                PLASMA_desc, descA,
                PLASMA_desc, descC,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        } else {
            plasma_dynamic_call_7(plasma_pzunmlq,
                PLASMA_enum, side,
                PLASMA_enum, trans,
                PLASMA_desc, descA,
                PLASMA_desc, descC,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
    }
    else {
        plasma_dynamic_call_8(plasma_pzunmlqrh,
            PLASMA_enum, side,
            PLASMA_enum, trans,
            PLASMA_desc, descA,
            PLASMA_desc, descC,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    return PLASMA_SUCCESS;
}
