/**
 *
 * @file zgetrs.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
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
 *  PLASMA_zgetrs - Solves a system of linear equations A * X = B, with a general N-by-N matrix A
 *  using the tile LU factorization computed by PLASMA_zgetrf.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = PlasmaNoTrans:   A * X = B     (No transpose)
 *          = PlasmaTrans:     A**T * X = B  (Transpose)
 *          = PlasmaConjTrans: A**H * X = B  (Conjugate transpose)
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_zgetrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_zgetrf.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \return <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrs_Tile
 * @sa PLASMA_zgetrs_Tile_Async
 * @sa PLASMA_cgetrs
 * @sa PLASMA_dgetrs
 * @sa PLASMA_sgetrs
 * @sa PLASMA_zgetrf
 *
 ******************************************************************************/
int PLASMA_zgetrs(PLASMA_enum trans, int N, int NRHS,
                  PLASMA_Complex64_t *A, int LDA,
                  int *IPIV,
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
        plasma_fatal_error("PLASMA_zgetrs", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( (trans != PlasmaNoTrans) &&
         (trans != PlasmaTrans)   &&
         (trans != PlasmaConjTrans)) {
        plasma_error("PLASMA_zgetrs", "illegal value of trans");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgetrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_zgetrs", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zgetrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_zgetrs", "illegal value of LDB");
        return -7;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGESV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrs", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,   
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zgetrs_Tile_Async(trans, &descA, IPIV, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
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
 *  PLASMA_zgetrs_Tile - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Tile equivalent of PLASMA_zgetrs().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended to specify the the form of the system of equations:
 *          = PlasmaNoTrans:   A * X = B     (No transpose)
 *          = PlasmaTrans:     A**T * X = B  (Transpose)
 *          = PlasmaConjTrans: A**H * X = B  (Conjugate transpose)
 *
 * @param[in] A
 *          The tile factors L and U from the factorization, computed by PLASMA_zgetrf.
 *
 * @param[in] IPIV
 *          The pivot indices from PLASMA_zgetrf.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrs
 * @sa PLASMA_zgetrs_Tile_Async
 * @sa PLASMA_cgetrs_Tile
 * @sa PLASMA_dgetrs_Tile
 * @sa PLASMA_sgetrs_Tile
 * @sa PLASMA_zgetrf_Tile
 *
 ******************************************************************************/
int PLASMA_zgetrs_Tile(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgetrs_Tile_Async(trans, A, IPIV, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgetrs_Tile_Async - Solves a system of linear equations using previously
 *  computed LU factorization.
 *  Non-blocking equivalent of PLASMA_zgetrs_Tile().
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
 * @sa PLASMA_zgetrs
 * @sa PLASMA_zgetrs_Tile
 * @sa PLASMA_cgetrs_Tile_Async
 * @sa PLASMA_dgetrs_Tile_Async
 * @sa PLASMA_sgetrs_Tile_Async
 * @sa PLASMA_zgetrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgetrs_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, int *IPIV, PLASMA_desc *B,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrs_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgetrs_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgetrs_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrs_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrs_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_zgetrs_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/

    if ( trans == PlasmaNoTrans )
    {
        plasma_dynamic_call_3(
            plasma_pzbarrier_tl2pnl,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*,  request);

        /* swap */
        plasma_dynamic_call_5(
            plasma_pzlaswp,
            PLASMA_desc, descB,
            int *,       IPIV,
            int,         1,
            PLASMA_sequence*, sequence,
            PLASMA_request*,  request);

        plasma_parallel_call_9(
            plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaLower,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        plasma_parallel_call_9(
            plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaUpper,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_parallel_call_9(
            plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaUpper,
            PLASMA_enum, trans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        plasma_parallel_call_9(
            plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaLower,
            PLASMA_enum, trans,
            PLASMA_enum, PlasmaUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descA,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        plasma_dynamic_call_3(
            plasma_pzbarrier_tl2pnl,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*,  request);

        /* swap */
        plasma_dynamic_call_5(
            plasma_pzlaswp,
            PLASMA_desc, descB,
            int *,       IPIV,
            int,         -1,
            PLASMA_sequence*, sequence,
            PLASMA_request*,  request);

        plasma_dynamic_call_3(
            plasma_pzbarrier_pnl2tl,
            PLASMA_desc, descB,
            PLASMA_sequence*, sequence,
            PLASMA_request*,  request);
    }
    return PLASMA_SUCCESS;
}
