/**
 *
 * @file zsyrk.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
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
 *  PLASMA_zsyrk - Performs one of the hermitian rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans  :   A is transposed.
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ) if trans == PlasmaNoTrans, otherwise LDA must
 *          be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zsyrk_Tile
 * @sa PLASMA_csyrk
 * @sa PLASMA_dsyrk
 * @sa PLASMA_ssyrk
 *
 ******************************************************************************/
int PLASMA_zsyrk(PLASMA_enum uplo, PLASMA_enum trans, int N, int K,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int LDC)
{
    int NB;
    int Am, An;
    int status;
    PLASMA_desc descA, descC;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zsyrk", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_zsyrk", "illegal value of uplo");
        return -1;
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_zsyrk", "illegal value of trans");
        return -2;
    }
    if ( trans == PlasmaNoTrans ) { 
        Am = N; An = K;
    } else {
        Am = K; An = N;
    }
    if (N < 0) {
        plasma_error("PLASMA_zsyrk", "illegal value of N");
        return -3;
    }
    if (K < 0) {
        plasma_error("PLASMA_zsyrk", "illegal value of K");
        return -4;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_zsyrk", "illegal value of LDA");
        return -7;
    }
    if (LDC < max(1, N)) {
        plasma_error("PLASMA_zsyrk", "illegal value of LDC");
        return -10;
    }

    /* Quick return */
    if (N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0.0) && beta == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZSYRK, N, K, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zsyrk", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N,   sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, 
                            sequence, &request);
        plasma_ziplap2tile( descC, C, NB, NB, LDC, N,  0, 0, N,  N, 
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zsyrk_Tile_Async(uplo, trans, alpha, &descA, beta, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, An,  sequence, &request);
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
 *  PLASMA_zsyrk_Tile - Performs rank k update.
 *  Tile equivalent of PLASMA_zsyrk().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of C is stored;
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = PlasmaNoTrans,
 *          and is N otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zsyrk_Tile
 * @sa PLASMA_csyrk
 * @sa PLASMA_dsyrk
 * @sa PLASMA_ssyrk
 *
 ******************************************************************************/
int PLASMA_zsyrk_Tile(PLASMA_enum uplo, PLASMA_enum trans,
                      PLASMA_Complex64_t alpha, PLASMA_desc *A,
                      PLASMA_Complex64_t beta,  PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zsyrk_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zsyrk_Tile_Async(uplo, trans, alpha, A, beta, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zsyrk_Tile_Async - Performs rank-k update.
 *  Non-blocking equivalent of PLASMA_zsyrk_Tile().
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
 * @sa PLASMA_zsyrk
 * @sa PLASMA_zsyrk_Tile
 * @sa PLASMA_csyrk_Tile_Async
 * @sa PLASMA_dsyrk_Tile_Async
 * @sa PLASMA_ssyrk_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zsyrk_Tile_Async(PLASMA_enum uplo, PLASMA_enum trans,
                            PLASMA_Complex64_t alpha, PLASMA_desc *A,
                            PLASMA_Complex64_t beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descC;
    int N, K;
    int Am, An, Amb;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zsyrk_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zsyrk_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zsyrk_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descC = *C;
    }
    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        plasma_error("PLASMA_zsyrk", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((trans != PlasmaNoTrans) && (trans != PlasmaTrans)) {
        plasma_error("PLASMA_zsyrk", "illegal value of transA");
        return plasma_request_fail(sequence, request, -2);
    }

    if ( trans == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
    }

    if (descC.mb != descC.nb) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Amb != descC.mb) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descC.m != descC.n) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "only square matrix C is supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (Am != descC.m) {
        plasma_error("PLASMA_zsyrk_Tile_Async", "sizes of matrices have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    N = descC.m;
    K = An;

    /* Quick return */
    if ( N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0) && beta == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    plasma_parallel_call_8(plasma_pzsyrk,
        PLASMA_enum, uplo,
        PLASMA_enum, trans,
        PLASMA_Complex64_t, alpha,
        PLASMA_desc, descA,
        PLASMA_Complex64_t, beta,
        PLASMA_desc, descC,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
