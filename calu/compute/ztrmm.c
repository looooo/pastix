/**
 *
 * @file ztrmm.c
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
 *  PLASMA_ztrmm - Computes B = alpha*op( A )*B or B = alpha*B*op( A ).
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B. NRHS >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_ztrmm_Tile
 * @sa PLASMA_ztrmm_Tile_Async
 * @sa PLASMA_ctrmm
 * @sa PLASMA_dtrmm
 * @sa PLASMA_strmm
 *
 ******************************************************************************/
int PLASMA_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                 PLASMA_enum transA, PLASMA_enum diag,
                 int N, int NRHS, PLASMA_Complex64_t alpha,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB)
{
    int NB, NA;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ztrmm", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if (side == PlasmaLeft) {
      NA = N;
    } else {
      NA = NRHS;
    }

    /* Check input arguments */
    if (side != PlasmaLeft && side != PlasmaRight) {
        plasma_error("PLASMA_ztrmm", "illegal value of side");
        return -1;
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_ztrmm", "illegal value of uplo");
        return -2;
    }
    if (transA != PlasmaConjTrans && 
        transA != PlasmaNoTrans   && 
        transA != PlasmaTrans ) 
    {
        plasma_error("PLASMA_ztrmm", "illegal value of transA");
        return -3;
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_ztrmm", "illegal value of diag");
        return -4;
    }
    if (N < 0) {
        plasma_error("PLASMA_ztrmm", "illegal value of N");
        return -5;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_ztrmm", "illegal value of NRHS");
        return -6;
    }
    if (LDA < max(1, NA)) {
        plasma_error("PLASMA_ztrmm", "illegal value of LDA");
        return -8;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_ztrmm", "illegal value of LDB");
        return -10;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZPOSV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ztrmm", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, NA,   0, 0, NA, NA,  
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N,  NRHS,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_ztrmm_Tile_Async(
        side, uplo, transA, diag, alpha, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, NA,    sequence, &request);
        plasma_ziptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
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
 *  PLASMA_ztrmm_Tile - Computes triangular solve.
 *  Tile equivalent of PLASMA_ztrmm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = PlasmaLeft:  A*X = B
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = PlasmaNoTrans:   A is transposed;
 *          = PlasmaTrans:     A is not transposed;
 *          = PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = PlasmaNonUnit: A is non unit;
 *          = PlasmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_ztrmm
 * @sa PLASMA_ztrmm_Tile_Async
 * @sa PLASMA_ctrmm_Tile
 * @sa PLASMA_dtrmm_Tile
 * @sa PLASMA_strmm_Tile
 *
 ******************************************************************************/
int PLASMA_ztrmm_Tile(PLASMA_enum side, PLASMA_enum uplo,
                      PLASMA_enum transA, PLASMA_enum diag,
                      PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ztrmm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_ztrmm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_ztrmm_Tile_Async - Performs triangular matrix multiplication.
 *  Non-blocking equivalent of PLASMA_ztrmm_Tile().
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
 * @sa PLASMA_ztrmm
 * @sa PLASMA_ztrmm_Tile
 * @sa PLASMA_ctrmm_Tile_Async
 * @sa PLASMA_dtrmm_Tile_Async
 * @sa PLASMA_strmm_Tile_Async
 *
 ******************************************************************************/
int PLASMA_ztrmm_Tile_Async(PLASMA_enum side, PLASMA_enum uplo,
                            PLASMA_enum transA, PLASMA_enum diag,
                            PLASMA_Complex64_t alpha, PLASMA_desc *A, PLASMA_desc *B,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_ztrmm_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_ztrmm_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_ztrmm_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ztrmm_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_ztrmm_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_ztrmm_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (side != PlasmaLeft && side != PlasmaRight) {
        plasma_error("PLASMA_ztrmm_Tile", "illegal value of side");
        return plasma_request_fail(sequence, request, -1);
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_ztrmm_Tile", "illegal value of uplo");
        return plasma_request_fail(sequence, request, -2);
    }
    if (transA != PlasmaConjTrans && transA != PlasmaNoTrans && transA != PlasmaTrans) {
        plasma_error("PLASMA_ztrmm_Tile", "illegal value of transA");
        return plasma_request_fail(sequence, request, -3);
    }
    if (diag != PlasmaUnit && diag != PlasmaNonUnit) {
        plasma_error("PLASMA_ztrmm_Tile", "illegal value of diag");
        return plasma_request_fail(sequence, request, -4);
    }

    /* Quick return */
    plasma_dynamic_call_9(plasma_pztrmm,
        PLASMA_enum, side,
        PLASMA_enum, uplo,
        PLASMA_enum, transA,
        PLASMA_enum, diag,
        PLASMA_Complex64_t, alpha,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
