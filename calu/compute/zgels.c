/**
 *
 * @file zgels.c
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
 *  PLASMA_zgels - solves overdetermined or underdetermined linear systems involving an M-by-N
 *  matrix A using the QR or the LQ factorization of A.  It is assumed that A has full rank.
 *  The following options are provided:
 *
 *  # trans = PlasmaNoTrans and M >= N: find the least squares solution of an overdetermined
 *    system, i.e., solve the least squares problem: minimize || B - A*X ||.
 *
 *  # trans = PlasmaNoTrans and M < N:  find the minimum norm solution of an underdetermined
 *    system A * X = B.
 *
 *  Several right hand side vectors B and solution vectors X can be handled in a single call;
 *  they are stored as the columns of the M-by-NRHS right hand side matrix B and the N-by-NRHS
 *  solution matrix X.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   the linear system involves A;
 *          = PlasmaConjTrans: the linear system involves A**H.
 *          Currently only PlasmaNoTrans is supported.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrices B and X.
 *          NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if M >= N, A is overwritten by details of its QR factorization as returned by
 *                     PLASMA_zgeqrf;
 *          if M < N, A is overwritten by details of its LQ factorization as returned by
 *                      PLASMA_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] descT
 *          On exit, auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS matrix B of right hand side vectors, stored columnwise;
 *          On exit, if return value = 0, B is overwritten by the solution vectors, stored
 *          columnwise:
 *          if M >= N, rows 1 to N of B contain the least squares solution vectors; the residual
 *          sum of squares for the solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column;
 *          if M < N, rows 1 to N of B contain the minimum norm solution vectors;
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= MAX(1,M,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgels_Tile
 * @sa PLASMA_zgels_Tile_Async
 * @sa PLASMA_cgels
 * @sa PLASMA_dgels
 * @sa PLASMA_sgels
 *
 ******************************************************************************/
int PLASMA_zgels(PLASMA_enum trans, int M, int N, int NRHS,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_desc *descT,
                 PLASMA_Complex64_t *B, int LDB)
{
    int i, j;
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgels", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (trans != PlasmaNoTrans) {
        plasma_error("PLASMA_zgels", "only PlasmaNoTrans supported");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (M < 0) {
        plasma_error("PLASMA_zgels", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgels", "illegal value of N");
        return -3;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_zgels", "illegal value of NRHS");
        return -4;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgels", "illegal value of LDA");
        return -6;
    }
    if (LDB < max(1, max(M, N))) {
        plasma_error("PLASMA_zgels", "illegal value of LDB");
        return -9;
    }
    /* Quick return */
    if (min(M, min(N, NRHS)) == 0) {
        for (i = 0; i < max(M, N); i++)
            for (j = 0; j < NRHS; j++)
                B[j*LDB+i] = 0.0;
        return PLASMA_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGELS, M, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgels", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( M >= N ) {
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            plasma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N, sequence, &request,
                                 plasma_desc_mat_free(&(descA)) );
            plasma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS, sequence, &request,
                                 plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
        } else {
            plasma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,
                                sequence, &request);
            plasma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS,
                                sequence, &request);
        }
    } else {
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            plasma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N, sequence, &request,
                                 plasma_desc_mat_free(&(descA)) );
            plasma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                                 plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)));
        } else {
            plasma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,
                                sequence, &request);
            plasma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,
                                sequence, &request);
        }
    }

    /* Call the tile interface */
    PLASMA_zgels_Tile_Async(PlasmaNoTrans, &descA, descT, &descB, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
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
 *  PLASMA_zgels_Tile - Solves overdetermined or underdetermined linear system of equations
 *  using the tile QR or the tile LQ factorization.
 *  Tile equivalent of PLASMA_zgels().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   the linear system involves A;
 *          = PlasmaConjTrans: the linear system involves A**H.
 *          Currently only PlasmaNoTrans is supported.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if M >= N, A is overwritten by details of its QR factorization as returned by
 *                     PLASMA_zgeqrf;
 *          if M < N, A is overwritten by details of its LQ factorization as returned by
 *                      PLASMA_zgelqf.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS matrix B of right hand side vectors, stored columnwise;
 *          On exit, if return value = 0, B is overwritten by the solution vectors, stored
 *          columnwise:
 *          if M >= N, rows 1 to N of B contain the least squares solution vectors; the residual
 *          sum of squares for the solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column;
 *          if M < N, rows 1 to N of B contain the minimum norm solution vectors;
 *
 *******************************************************************************
 *
 * @return
 *          \return PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgels
 * @sa PLASMA_zgels_Tile_Async
 * @sa PLASMA_cgels_Tile
 * @sa PLASMA_dgels_Tile
 * @sa PLASMA_sgels_Tile
 *
 ******************************************************************************/
int PLASMA_zgels_Tile(PLASMA_enum trans, PLASMA_desc *A,
                      PLASMA_desc *T, PLASMA_desc *B)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgels_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgels_Tile_Async(trans, A, T, B, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgels_Tile_Async - Solves overdetermined or underdetermined linear
 *  system of equations using the tile QR or the tile LQ factorization.
 *  Non-blocking equivalent of PLASMA_zgels_Tile().
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
 * @sa PLASMA_zgels
 * @sa PLASMA_zgels_Tile
 * @sa PLASMA_cgels_Tile_Async
 * @sa PLASMA_dgels_Tile_Async
 * @sa PLASMA_sgels_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgels_Tile_Async(PLASMA_enum trans, PLASMA_desc *A,
                            PLASMA_desc *T, PLASMA_desc *B,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descB;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgels_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgels_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgels_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgels_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgels_Tile", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgels_Tile", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb) {
        plasma_error("PLASMA_zgels_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (trans != PlasmaNoTrans) {
        plasma_error("PLASMA_zgels_Tile", "only PlasmaNoTrans supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
    }
    /* Quick return  - currently NOT equivalent to LAPACK's:
    if (min(M, min(N, NRHS)) == 0) {
        for (i = 0; i < max(M, N); i++)
            for (j = 0; j < NRHS; j++)
                B[j*LDB+i] = 0.0;
        return PLASMA_SUCCESS;
    }
*/
    if (descA.m >= descA.n) {
        if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
            plasma_parallel_call_4(plasma_pzgeqrf,
                PLASMA_desc, descA,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);

            plasma_parallel_call_7(plasma_pzunmqr,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        else {
            plasma_dynamic_call_5(plasma_pzgeqrfrh,
                PLASMA_desc, descA,
                PLASMA_desc, descT,
                PLASMA_enum, PLASMA_RHBLK,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);

            plasma_dynamic_call_8(plasma_pzunmqrrh,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_enum, PLASMA_RHBLK,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        plasma_parallel_call_9(plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaUpper,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, plasma_desc_submatrix(descA, 0, 0, descA.n, descA.n),
            PLASMA_desc, plasma_desc_submatrix(descB, 0, 0, descA.n, descB.n),
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);
    }
    else {
        plasma_parallel_call_3(plasma_pztile_zero,
            PLASMA_desc, plasma_desc_submatrix(descB, descA.m, 0, descA.n-descA.m, descB.n),
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
            plasma_parallel_call_4(plasma_pzgelqf,
                PLASMA_desc, descA,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        else {
            plasma_dynamic_call_5(plasma_pzgelqfrh,
                PLASMA_desc, descA,
                PLASMA_desc, descT,
                PLASMA_enum, PLASMA_RHBLK,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        plasma_parallel_call_9(plasma_pztrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaLower,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, plasma_desc_submatrix(descA, 0, 0, descA.m, descA.m),
            PLASMA_desc, plasma_desc_submatrix(descB, 0, 0, descA.m, descB.n),
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) {
            plasma_parallel_call_7(plasma_pzunmlq,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
        else {
            plasma_dynamic_call_8(plasma_pzunmlqrh,
                PLASMA_enum, PlasmaLeft,
                PLASMA_enum, PlasmaConjTrans,
                PLASMA_desc, descA,
                PLASMA_desc, descB,
                PLASMA_desc, descT,
                PLASMA_enum, PLASMA_RHBLK,
                PLASMA_sequence*, sequence,
                PLASMA_request*, request);
        }
    }
    return PLASMA_SUCCESS;
}
