/**
 *
 * @file zungqr.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zungqr - Generates an M-by-N matrix Q with orthonormal columns, which is defined as the
 *  first N columns of a product of the elementary reflectors returned by PLASMA_zgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix Q. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix Q. N >= M.
 *
 * @param[in] K
 *          The number of columns of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_zgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by PLASMA_zgeqrf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zungqr_Tile
 * @sa PLASMA_zungqr_Tile_Async
 * @sa PLASMA_cungqr
 * @sa PLASMA_dorgqr
 * @sa PLASMA_sorgqr
 * @sa PLASMA_zgeqrf
 *
 ******************************************************************************/
int PLASMA_zungqr(int M, int N, int K,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_desc *descT,
                  PLASMA_Complex64_t *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descQ;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zungqr", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zungqr", "illegal value of M");
        return -1;
    }
    if (N < 0 || N > M) {
        plasma_error("PLASMA_zungqr", "illegal value of N");
        return -2;
    }
    if (K < 0 || K > N) {
        plasma_error("PLASMA_zungqr", "illegal value of K");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zungqr", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, M)) {
        plasma_error("PLASMA_zungqr", "illegal value of LDQ");
        return -8;
    }
    if (min(M, min(N, K)) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGELS, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zungqr", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, K, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descQ)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, K,
                            sequence, &request);
        plasma_ziplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zungqr_Tile_Async(&descA, descT, &descQ, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descQ);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, K,  sequence, &request);
        plasma_ziptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);
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
 *  PLASMA_zungqr_Tile - Generates an M-by-N matrix Q with orthonormal columns, which is defined as the
 *  first N columns of a product of the elementary reflectors returned by PLASMA_zgeqrf.
 *  All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the QR factorization of the original matrix A as returned by PLASMA_zgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by PLASMA_zgeqrf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zungqr
 * @sa PLASMA_zungqr_Tile_Async
 * @sa PLASMA_cungqr_Tile
 * @sa PLASMA_dorgqr_Tile
 * @sa PLASMA_sorgqr_Tile
 * @sa PLASMA_zgeqrf_Tile
 *
 ******************************************************************************/
int PLASMA_zungqr_Tile(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zungqr_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zungqr_Tile_Async(A, T, Q, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  Non-blocking equivalent of PLASMA_zungqr_Tile().
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
 * @sa PLASMA_zungqr
 * @sa PLASMA_zungqr_Tile
 * @sa PLASMA_cungqr_Tile_Async
 * @sa PLASMA_dorgqr_Tile_Async
 * @sa PLASMA_sorgqr_Tile_Async
 * @sa PLASMA_zgeqrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zungqr_Tile_Async(PLASMA_desc *A, PLASMA_desc *T, PLASMA_desc *Q,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descQ;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zungqr_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zungqr_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zungqr_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zungqr_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zungqr_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    if (plasma_desc_check(Q) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zungqr_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descQ = *Q;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descQ.nb != descQ.mb) {
        plasma_error("PLASMA_zungqr_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (N <= 0)
        return PLASMA_SUCCESS;
*/
    if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
        plasma_dynamic_call_5(plasma_pzungqr,
            PLASMA_desc, descA,
            PLASMA_desc, descQ,
            PLASMA_desc, descT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_dynamic_call_6(plasma_pzungqrrh,
            PLASMA_desc, descA,
            PLASMA_desc, descQ,
            PLASMA_desc, descT,
            PLASMA_enum, PLASMA_RHBLK,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }

    return PLASMA_SUCCESS;
}
