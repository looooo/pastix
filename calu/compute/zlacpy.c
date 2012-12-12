/**
 *
 * @file zlacpy.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zlacpy copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part
 *            = PlasmaLower: Lower triangular part
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0. 
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = PlasmaUpper, only the upper trapezium
 *          is accessed; if UPLO = PlasmaLower, only the lower trapezium is
 *          accessed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa PLASMA_zlacpy_Tile
 * @sa PLASMA_zlacpy_Tile_Async
 * @sa PLASMA_clacpy
 * @sa PLASMA_dlacpy
 * @sa PLASMA_slacpy
 *
 ******************************************************************************/
int PLASMA_zlacpy(PLASMA_enum uplo, int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_Complex64_t *B, int LDB)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlacpy", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (uplo != PlasmaUpperLower) && 
         (uplo != PlasmaUpper) &&
         (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_zlacpy", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        plasma_error("PLASMA_zlacpy", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_zlacpy", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zlacpy", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, M)) {
        plasma_error("PLASMA_zlacpy", "illegal value of LDB");
        return -7;
    }

    /* Quick return */
    if (min(N, M) == 0)
      return (double)0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zlacpy", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)) );
    } else {
        plasma_ziplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
        plasma_ziplap2tile(  descB, B, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zlacpy_Tile_Async(uplo, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ziptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
        plasma_dynamic_sync();
    }

    plasma_sequence_destroy(plasma, sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zlacpy_Tile - Tile equivalent of PLASMA_zlacpy().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = PlasmaUpperLower: All the matrix A
 *            = PlasmaUpper: Upper triangular part
 *            = PlasmaLower: Lower triangular part
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = PlasmaUpper, only the upper trapezium
 *          is accessed; if UPLO = PlasmaLower, only the lower trapezium is
 *          accessed.
 *
 * @param[out] B
 *          The M-by-N matrix B.
 *          On exit, B = A in the locations specified by UPLO.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zlacpy
 * @sa PLASMA_zlacpy_Tile_Async
 * @sa PLASMA_clacpy_Tile
 * @sa PLASMA_dlacpy_Tile
 * @sa PLASMA_slacpy_Tile
 *
 ******************************************************************************/
int PLASMA_zlacpy_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlacpy_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zlacpy_Tile_Async(uplo, A, B, sequence, &request);
    plasma_dynamic_sync();
    plasma_sequence_destroy(plasma, sequence);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zlacpy_Tile_Async - Non-blocking equivalent of PLASMA_zlacpy_Tile().
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
 * @sa PLASMA_zlacpy
 * @sa PLASMA_zlacpy_Tile
 * @sa PLASMA_clacpy_Tile_Async
 * @sa PLASMA_dlacpy_Tile_Async
 * @sa PLASMA_slacpy_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zlacpy_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zlacpy_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zlacpy_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zlacpy_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zlacpy_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zlacpy_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zlacpy_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( (uplo != PlasmaUpperLower) && 
         (uplo != PlasmaUpper) &&
         (uplo != PlasmaLower) ) {
        plasma_error("PLASMA_zlacpy_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (min(descA.m, descA.n) == 0) {
        return PLASMA_SUCCESS;
    }

    plasma_parallel_call_5(plasma_pzlacpy,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
