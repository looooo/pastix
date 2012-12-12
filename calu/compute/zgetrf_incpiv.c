/**
 *
 * @file zgetrf_incpiv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
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
 *  PLASMA_zgetrf_incpiv - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descL
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by PLASMA_zgetrs_incpiv to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrf_incpiv_Tile
 * @sa PLASMA_zgetrf_incpiv_Tile_Async
 * @sa PLASMA_cgetrf_incpiv
 * @sa PLASMA_dgetrf_incpiv
 * @sa PLASMA_sgetrf_incpiv
 * @sa PLASMA_zgetrs_incpiv
 *
 ******************************************************************************/
int PLASMA_zgetrf_incpiv(int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_desc *descL, int *IPIV)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_incpiv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zgetrf_incpiv", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgetrf_incpiv", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgetrf_incpiv", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGESV, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrf_incpiv", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zgetrf_incpiv_Tile_Async(&descA, descL, IPIV, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
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
 *  PLASMA_zgetrf_incpiv_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of PLASMA_zgetrf_incpiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] L
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          required by PLASMA_zgetrs_incpiv to solve the system of equations.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrf_incpiv
 * @sa PLASMA_zgetrf_incpiv_Tile_Async
 * @sa PLASMA_cgetrf_incpiv_Tile
 * @sa PLASMA_dgetrf_incpiv_Tile
 * @sa PLASMA_sgetrf_incpiv_Tile
 * @sa PLASMA_zgetrs_incpiv_Tile
 *
 ******************************************************************************/
int PLASMA_zgetrf_incpiv_Tile(PLASMA_desc *A, PLASMA_desc *L, int *IPIV)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_incpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgetrf_incpiv_Tile_Async(A, L, IPIV, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgetrf_incpiv_Tile_Async - Computes the tile LU factorization of a matrix.
 *  Non-blocking equivalent of PLASMA_zgetrf_incpiv_Tile().
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
 * @sa PLASMA_zgetrf_incpiv
 * @sa PLASMA_zgetrf_incpiv_Tile
 * @sa PLASMA_cgetrf_incpiv_Tile_Async
 * @sa PLASMA_dgetrf_incpiv_Tile_Async
 * @sa PLASMA_sgetrf_incpiv_Tile_Async
 * @sa PLASMA_zgetrs_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgetrf_incpiv_Tile_Async(PLASMA_desc *A, PLASMA_desc *L, int *IPIV,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descL;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_incpiv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_incpiv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_incpiv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrf_incpiv_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(L) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrf_incpiv_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descL = *L;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgetrf_incpiv_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;
*/

    plasma_parallel_call_5(plasma_pzgetrf_incpiv,
        PLASMA_desc, descA,
        PLASMA_desc, descL,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
