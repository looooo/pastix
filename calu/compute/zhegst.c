/**
 *
 * @file zhegst.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zhegst - reduces a complex Hermitian-definite generalized
 *  eigenproblem to standard form.
 *  If PlasmaItype == 1, the problem is A*x = lambda*B*x, and A is
 *  overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
 *  If PlasmaItype == 2 or 3, the problem is A*B*x = lambda*x or B*A*x
 *  = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.  B must
 *  have been previously factorized as U**H*U or L*L**H by
 *  PLASMA_ZPOTRF.
 *
 *******************************************************************************
 *
 * @param[in] itype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x
 *          = 3: B*A*x=(lambda)*x
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The order of the matrices A and B. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value == 0, the transformed matrix,
 *          stored in the same format as A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the triangular factor from the Cholesky
 *          factorization of B, as returned by PLASMA_ZPOTRF.
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
 * @sa PLASMA_zhegst_Tile
 * @sa PLASMA_zhegst_Tile_Async
 * @sa PLASMA_chegst
 * @sa PLASMA_dsygst
 * @sa PLASMA_ssygst
 *
 ******************************************************************************/
int PLASMA_zhegst(PLASMA_enum itype, PLASMA_enum uplo, int N,
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
        plasma_fatal_error("PLASMA_zhegst", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (itype != 1 && itype != 2 && itype != 3) {
        plasma_error("PLASMA_zhegst", "Illegal value of itype");
        return -1;
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_zhegst", "Illegal value of uplo");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_zhegst", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zhegst", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_zhegst", "illegal value of LDB");
        return -7;
    }
    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZHEGST, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegst", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descB)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zhegst_Tile_Async(itype, uplo, &descA, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_zooptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_ziptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);
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
 *  PLASMA_zhegst_Tile - reduces a complex Hermitian-definite
 *  generalized eigenproblem to standard form.
 *  If PlasmaItype == 1, the problem is A*x = lambda*B*x, and A is
 *  overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
 *  If PlasmaItype == 2 or 3, the problem is A*B*x = lambda*x or B*A*x
 *  = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.  B must
 *  have been previously factorized as U**H*U or L*L**H by
 *  PLASMA_ZPOTRF.
 *  ONLY PlasmaItype == 1 and PlasmaLower supported!
 *  Tile equivalent of PLASMA_zhegst().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] itype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x
 *          = 3: B*A*x=(lambda)*x
 *          Currently only PlasmaItype == 1 is supported.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *          Currently only PlasmaLower is supported.
 *
 * @param[in,out] A
 *          On entry, the symmetric (or Hermitian) matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of A contains the upper triangular part of the matrix
 *          A, and the strictly lower triangular part of A is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of A contains the lower triangular part of the matrix
 *          A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value == 0, the transformed matrix,
 *          stored in the same format as A.
 *
 * @param[in,out] B
 *          On entry, the triangular factor from the Cholesky
 *          factorization of B, as returned by PLASMA_ZPOTRF.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zhegst
 * @sa PLASMA_zhegst_Tile_Async
 * @sa PLASMA_chegst_Tile
 * @sa PLASMA_dsygst_Tile
 * @sa PLASMA_ssygst_Tile
 * @sa PLASMA_zhegst_Tile
 *
 ******************************************************************************/
int PLASMA_zhegst_Tile(PLASMA_enum itype, PLASMA_enum uplo,
                       PLASMA_desc *A,
                       PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhegst_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zhegst_Tile_Async(itype, uplo, A, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zhegst_Tile_Async - reduces a complex Hermitian-definite
 *  generalized eigenproblem to standard form.
 *  If PlasmaItype == 1, the problem is A*x = lambda*B*x, and A is
 *  overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
 *  If PlasmaItype == 2 or 3, the problem is A*B*x = lambda*x or B*A*x
 *  = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.  B must
 *  have been previously factorized as U**H*U or L*L**H by
 *  PLASMA_ZPOTRF.
 *  ONLY PlasmaItype == 1 and PlasmaLower supported!
 *  Non-blocking equivalent of PLASMA_zhegst_Tile().
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
 * @sa PLASMA_zhegst
 * @sa PLASMA_zhegst_Tile
 * @sa PLASMA_chegst_Tile_Async
 * @sa PLASMA_dsygst_Tile_Async
 * @sa PLASMA_ssygst_Tile_Async
 * @sa PLASMA_zhegv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zhegst_Tile_Async(PLASMA_enum itype, PLASMA_enum uplo,
                             PLASMA_desc *A,
                             PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhegst_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zhegst_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zhegst_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegst_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegst_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zhegst_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }


    /*
     * Transform Hermitian-definite generalized eigenproblem
     * to standard form
     */
    plasma_dynamic_call_6(plasma_pzhegst,
        PLASMA_enum, itype,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_desc, descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
