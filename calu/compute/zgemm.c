/**
 *
 * @file zgemm.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Emmanuel Agullo
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zgemm - Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = PlasmaNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgemm_Tile
 * @sa PLASMA_cgemm
 * @sa PLASMA_dgemm
 * @sa PLASMA_sgemm
 *
 ******************************************************************************/
int PLASMA_zgemm(PLASMA_enum transA, PLASMA_enum transB, int M, int N, int K,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA,
                                           PLASMA_Complex64_t *B, int LDB,
                 PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int LDC)
{
    int NB;
    int Am, An, Bm, Bn;
    int status;
    PLASMA_desc descA, descB, descC;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgemm", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        plasma_error("PLASMA_zgemm", "illegal value of transA");
        return -1;
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        plasma_error("PLASMA_zgemm", "illegal value of transB");
        return -2;
    }
    if ( transA == PlasmaNoTrans ) { 
        Am = M; An = K;
    } else {
        Am = K; An = M;
    }
    if ( transB == PlasmaNoTrans ) { 
        Bm = K; Bn = N;
    } else {
        Bm = N; Bn = K;
    }
    if (M < 0) {
        plasma_error("PLASMA_zgemm", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgemm", "illegal value of N");
        return -4;
    }
    if (K < 0) {
        plasma_error("PLASMA_zgemm", "illegal value of N");
        return -5;
    }
    if (LDA < max(1, Am)) {
        plasma_error("PLASMA_zgemm", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, Bm)) {
        plasma_error("PLASMA_zgemm", "illegal value of LDB");
        return -10;
    }
    if (LDC < max(1, M)) {
        plasma_error("PLASMA_zgemm", "illegal value of LDC");
        return -13;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0) && beta == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm", "plasma_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
        plasma_zooplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N,  sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descC)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, An, 0, 0, Am, An, 
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, Bn, 0, 0, Bm, Bn, 
                            sequence, &request);
        plasma_ziplap2tile( descC, C, NB, NB, LDC, N,  0, 0, M,  N,  
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zgemm_Tile_Async(
        transA, transB, alpha, &descA, &descB, beta, &descC, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descC, C, NB, NB, LDC, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
        plasma_desc_mat_free(&descC);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, An,  sequence, &request);
        plasma_ziptile2lap( descB, B, NB, NB, LDB, Bn,  sequence, &request);
        plasma_ziptile2lap( descC, C, NB, NB, LDC, N,   sequence, &request);
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
 *  PLASMA_zgemm_Tile - Performs matrix multiplication.
 *  Tile equivalent of PLASMA_zgemm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = PlasmaNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgemm
 * @sa PLASMA_zgemm_Tile_Async
 * @sa PLASMA_cgemm_Tile
 * @sa PLASMA_dgemm_Tile
 * @sa PLASMA_sgemm_Tile
 *
 ******************************************************************************/
int PLASMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB,
                      PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B,
                      PLASMA_Complex64_t beta,  PLASMA_desc *C)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgemm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgemm_Tile_Async(transA, transB, alpha, A, B, beta, C, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of PLASMA_zgemm_Tile().
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
 * @sa PLASMA_zgemm
 * @sa PLASMA_zgemm_Tile
 * @sa PLASMA_cgemm_Tile_Async
 * @sa PLASMA_dgemm_Tile_Async
 * @sa PLASMA_sgemm_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB,
                            PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B,
                            PLASMA_Complex64_t beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descC;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgemm_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgemm_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgemm_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    if (plasma_desc_check(C) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgemm_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descC = *C;
    }
    /* Check input arguments */
    if ((transA != PlasmaNoTrans) && (transA != PlasmaTrans) && (transA != PlasmaConjTrans)) {
        plasma_error("PLASMA_zgemm_Tile_Async", "illegal value of transA");
        return plasma_request_fail(sequence, request, -1);
    }
    if ((transB != PlasmaNoTrans) && (transB != PlasmaTrans) && (transB != PlasmaConjTrans)) {
        plasma_error("PLASMA_zgemm_Tile_Async", "illegal value of transB");
        return plasma_request_fail(sequence, request, -2);
    }

    if ( transA == PlasmaNoTrans ) {
        Am  = descA.m;
        An  = descA.n;
        Amb = descA.mb;
        Anb = descA.nb;
        Ai  = descA.i;
        Aj  = descA.j;
    } else {
        Am  = descA.n;
        An  = descA.m;
        Amb = descA.nb;
        Anb = descA.mb;
        Ai  = descA.j;
        Aj  = descA.i;
    }

    if ( transB == PlasmaNoTrans ) {
        Bm  = descB.m;
        Bn  = descB.n;
        Bmb = descB.mb;
        Bnb = descB.nb;
        Bi  = descB.i;
        Bj  = descB.j;
    } else {
        Bm  = descB.n;
        Bn  = descB.m;
        Bmb = descB.nb;
        Bnb = descB.mb;
        Bi  = descB.j;
        Bj  = descB.i;
    }

    if ( (Amb != descC.mb) || (Anb != Bmb) || (Bnb != descC.nb) ) {
        plasma_error("PLASMA_zgemm_Tile_Async", "tile sizes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Am != descC.m) || (An != Bm) || (Bn != descC.n) ) {
        plasma_error("PLASMA_zgemm_Tile_Async", "sizes of matrices have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (Ai != descC.i) || (Aj != Bi) || (Bj != descC.j) ) {
        plasma_error("PLASMA_zgemm_Tile_Async", "start indexes have to match");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    M = descC.m;
    N = descC.n;
    K = An;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == (PLASMA_Complex64_t)0.0 || K == 0) && beta == (PLASMA_Complex64_t)1.0))
        return PLASMA_SUCCESS;

    plasma_parallel_call_9(plasma_pzgemm,
        PLASMA_enum, transA,
        PLASMA_enum, transB,
        PLASMA_Complex64_t, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_Complex64_t, beta,
        PLASMA_desc, descC,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
