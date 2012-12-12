/**
 *
 * @file zgebrd.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Azzam Haidar
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
 *  PLASMA_zgebrd - reduces a general complex M-by-N matrix A to upper or lower
 *  bidiagonal form B using a two-stage approach
 *  First stage: reduction to band bidiagonal form (orthogonal matrices Q1 and P1);
 *  Second stage: reduction from band to bidiagonal form (orthogonal matrices
 *  Q2 and P2).
 *  Let Q = Q1 * Q2 be the global left unitary transformation;
 *  Let P = P1 * P2 be the global right unitary transformation;
 *  Q**H * A * P = B.
 *  If M >= N, B is upper bidiagonal; if M < N, B is lower bidiagonal.
 *  Not LAPACK Compliant for now!
 *  Note: T is incomplete and contains only the block reflectors of the first stage.
 *  Therefore, Q and P can not be built completely.
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
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if M >= N, the diagonal and the first superdiagonal are
 *            overwritten with the upper bidiagonal matrix B; the
 *            elements below the diagonal, with the array T, represent
 *            the unitary matrix Q as a product of elementary
 *            reflectors, and the elements above the first superdiagonal,
 *            with the array T, represent the unitary matrix P as
 *            a product of elementary reflectors;
 *          if M < N, the diagonal and the first subdiagonal are
 *            overwritten with the lower bidiagonal matrix B; the
 *            elements below the first subdiagonal, with the array T,
 *            represent the unitary matrix Q as a product of
 *            elementary reflectors, and the elements above the diagonal,
 *            with the array T, represent the unitary matrix P as
 *            a product of elementary reflectors.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] D
 *          On exit, the diagonal elements of the bidiagonal matrix:
 *          D(i) = A(i,i).
 *          Dimension (min(M,N)).
 *
 * @param[out] E
 *          On exit, the off-diagonal elements of the bidiagonal matrix:
 *          if M >= N, E(i) = A(i,i+1) for i = 1,2,...,N-1;
 *          if M < N, E(i) = A(i+1,i) for i = 1,2,...,M-1.
 *          Dimension (min(M,N)-1).
 *
 * @param[out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgebrd_Tile
 * @sa PLASMA_zgebrd_Tile_Async
 * @sa PLASMA_cgebrd
 * @sa PLASMA_dgebrd
 * @sa PLASMA_sgebrd
 *
 ******************************************************************************/
int PLASMA_zgebrd(int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  double *D,
                  double *E,
                  PLASMA_desc *descT)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgebrd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zgebrd", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgebrd", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgebrd", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (min(M, N) == 0) {
        return PLASMA_SUCCESS;
    }

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGEBRD, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgebrd", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N,
                            sequence, &request);
    }

    /* Call the tile interface */
    PLASMA_zgebrd_Tile_Async(&descA, D, E, descT, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB,  LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB,  LDA, N,  sequence, &request);
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
 *  PLASMA_zgebrd_Tile - reduces a general complex M-by-N matrix A to upper or lower
 *  bidiagonal form B using a two-stage approach
 *  First stage: reduction to band bidiagonal form (orthogonal matrices Q1 and P1);
 *  Second stage: reduction from band to bidiagonal form (orthogonal matrices
 *  Q2 and P2).
 *  Let Q = Q1 * Q2 be the global left unitary transformation;
 *  Let P = P1 * P2 be the global right unitary transformation;
 *  Q**H * A * P = B.
 *  If M >= N, B is upper bidiagonal; if M < N, B is lower bidiagonal.
 *  Note: T is incomplete and contains only the block reflectors of the first stage.
 *  Therefore, Q and P can not be built completely.
 *  Tile equivalent of PLASMA_zgebrd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if M >= N, the diagonal and the first superdiagonal are
 *            overwritten with the upper bidiagonal matrix B; the
 *            elements below the diagonal, with the array T, represent
 *            the unitary matrix Q as a product of elementary
 *            reflectors, and the elements above the first superdiagonal,
 *            with the array T, represent the unitary matrix P as
 *            a product of elementary reflectors;
 *          if M < N, the diagonal and the first subdiagonal are
 *            overwritten with the lower bidiagonal matrix B; the
 *            elements below the first subdiagonal, with the array T,
 *            represent the unitary matrix Q as a product of
 *            elementary reflectors, and the elements above the diagonal,
 *            with the array T, represent the unitary matrix P as
 *            a product of elementary reflectors.
 *
 * @param[out] D
 *          The double precision array containing the diagonal elements
 *          of the bidiagonal matrix B:
 *          D(i) = A(i,i).
 *          Dimension (min(M,N)).
 *
 * @param[out] E
 *          The double precision array containing the off-diagonal elements
 *          of the bidiagonal matrix B:
 *          if M >= N, E(i) = A(i,i+1) for i = 1,2,...,N-1;
 *          if M < N, E(i) = A(i+1,i) for i = 1,2,...,M-1.
 *          Dimension (min(M,N)-1).
 *
 * @param[out] T
 *          On exit, contains auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \return PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgebrd
 * @sa PLASMA_zgebrd_Tile_Async
 * @sa PLASMA_cgebrd_Tile
 * @sa PLASMA_dgebrd_Tile
 * @sa PLASMA_sgebrd_Tile
 *
 ******************************************************************************/
int PLASMA_zgebrd_Tile(PLASMA_desc *A,
                      double *D, double *E, PLASMA_desc *T)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgebrd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgebrd_Tile_Async(A, D, E, T, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgebrd_Tile_Async - reduces a general complex M-by-N matrix A to upper or lower
 *  bidiagonal form B using a two-stage approach
 *  First stage: reduction to band bidiagonal form (orthogonal matrices Q1 and P1);
 *  Second stage: reduction from band to bidiagonal form (orthogonal matrices
 *  Q2 and P2).
 *  Let Q = Q1 * Q2 be the global left unitary transformation;
 *  Let P = P1 * P2 be the global right unitary transformation;
 *  Q**H * A * P = B.
 *  If M >= N, B is upper bidiagonal; if M < N, B is lower bidiagonal.
 *  Note: T is incomplete and contains only the block reflectors of the first stage.
 *  Therefore, Q and P can not be built completely.
 *  Non-blocking equivalent of PLASMA_zgebrd_Tile().
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
 * @sa PLASMA_zgebrd
 * @sa PLASMA_zgebrd_Tile
 * @sa PLASMA_cgebrd_Tile_Async
 * @sa PLASMA_dgebrd_Tile_Async
 * @sa PLASMA_sgebrd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgebrd_Tile_Async(PLASMA_desc *A,
                             double *D, double *E, PLASMA_desc *T,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA ;
    PLASMA_desc descT ;

    plasma_context_t *plasma;
    plasma = plasma_context_self();

    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgebrd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgebrd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgebrd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgebrd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgebrd_Tile_Async", "invalid fifth descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgebrd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /*
     * Reduction to BAND bidiagonal form
     * May be further optimized using the algo described in Trefethen
     */
    /* if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
    plasma_dynamic_call_4(plasma_pzgebrd_ge2tb,
                          PLASMA_desc, descA,
                          PLASMA_desc, descT,
                          PLASMA_sequence*, sequence,
                          PLASMA_request*, request);
    /* } */
    /* else { */
    /*     plasma_dynamic_call_4(plasma_pzgebrd_ge2tb_rh, */
    /*         PLASMA_desc, descA, */
    /*         PLASMA_desc, descT, */
    /*         PLASMA_sequence*, sequence, */
    /*         PLASMA_request*, request); */
    /* } */

    /*
     * Set the V's to zero before the 2nd stage i.e., bulge chasing
     */
    plasma_dynamic_call_5(plasma_pzlaset2,
        PLASMA_enum, PlasmaLower,
        PLASMA_Complex64_t, 0.0,
        PLASMA_desc, descA.m >= descA.n ? descA : plasma_desc_submatrix(descA, descA.mb, 0, descA.m-descA.mb, descA.n),
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_5(plasma_pzlaset2,
        PLASMA_enum, PlasmaUpper,
        PLASMA_Complex64_t, 0.0,
        PLASMA_desc, descA.m >= descA.n ? plasma_desc_submatrix(descA, 0, descA.nb, descA.m, descA.n-descA.nb) : descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /*
     * Reduction from BAND bidiagonal to the final condensed form
     */
    plasma_dynamic_call_7(plasma_pzgebrd_tb2bd,
        PLASMA_enum, descA.m >= descA.n ? PlasmaUpper : PlasmaLower,
        PLASMA_desc, descA,
        double*, D,
        double*, E,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
