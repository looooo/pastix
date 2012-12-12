/**
 *
 * @file zgeqp3.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zgeqp3 - Computes the QR factorization with column pivoting
 *                  of a complex M-by-N matrix A: A*P = Q*R.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] jpvt
 *          Integer array of dimension N.
 *          On exit, if jpvt(j)=k, then the j-th column of A*P was the k-th column of A.
 *
 * @param[out] tau
 *          On exit, scalars that define Householder reflectors, size n.
 *
 * @param[out] work
 *          Workspace of size (n + 1)*nb.
 *
 * @param[out] rwork
 *          Workspace of size 2*n.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, a numerical instability was detected at column i.
 *                     This case is not currently handled, but in the future
 *                     it will be handled internally.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgeqp3_Tile
 * @sa PLASMA_zgeqp3_Tile_Async
 * @sa PLASMA_cgeqp3
 * @sa PLASMA_dgeqp3
 * @sa PLASMA_sgeqp3
 *
 ******************************************************************************/
int PLASMA_zgeqp3( int M, int N,
                   PLASMA_Complex64_t *A, int LDA,
                   int *jpvt, PLASMA_Complex64_t *tau,
                   PLASMA_Complex64_t *work, double *rwork )
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zgeqp3", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgeqp3", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgeqp3", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGELS, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgeqp3", "plasma_tune() failed");
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
    PLASMA_zgeqp3_Tile_Async(&descA, jpvt, tau, work, rwork, sequence, &request);

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
 *  PLASMA_zgeqp3_Tile - Computes the tile QR factorization
 *                       with column pivoting of a matrix.
 *  Tile equivalent of PLASMA_zgeqp3().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[in,out] jpvt
 *          Integer array of dimension N.
 *          On exit, if jpvt(j)=k, then the j-th column of A*P was the k-th column of A.
 *
 * @param[out] tau
 *          On exit, scalars that define Householder reflectors, size n.
 *
 * @param[out] work
 *          Workspace of size (n + 1)*nb.
 *
 * @param[out] rwork
 *          Workspace of size 2*n.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, a numerical instability was detected at column i.
 *                     This case is not currently handled, but in the future
 *                     it will be handled internally.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgeqp3
 * @sa PLASMA_zgeqp3_Tile_Async
 * @sa PLASMA_cgeqp3_Tile
 * @sa PLASMA_dgeqp3_Tile
 * @sa PLASMA_sgeqp3_Tile
 *
 ******************************************************************************/
int PLASMA_zgeqp3_Tile( PLASMA_desc *A,
                        int *jpvt, PLASMA_Complex64_t *tau,
                        PLASMA_Complex64_t *work, double *rwork )
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgeqp3_Tile_Async(A, jpvt, tau, work, rwork, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgeqp3_Tile_Async - Computes the tile QR factorization
 *                             with column pivoting of a matrix.
 *  Non-blocking equivalent of PLASMA_zgeqp3_Tile().
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
 * @sa PLASMA_zgeqp3
 * @sa PLASMA_zgeqp3_Tile
 * @sa PLASMA_cgeqp3_Tile_Async
 * @sa PLASMA_dgeqp3_Tile_Async
 * @sa PLASMA_sgeqp3_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgeqp3_Tile_Async( PLASMA_desc *A, int *jpvt, PLASMA_Complex64_t *tau,
                              PLASMA_Complex64_t *work, double *rwork,
                              PLASMA_sequence *sequence, PLASMA_request *request )
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_zgeqp3_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgeqp3_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (jpvt == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL jpvt");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (tau == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL tau");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (work == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL work");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (rwork == NULL) {
        plasma_fatal_error("PLASMA_zgeqp3_Tile", "NULL rwork");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgeqp3_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    plasma_dynamic_call_7(plasma_pzgeqp3,
        PLASMA_desc,         descA,
        int*,                jpvt,
        PLASMA_Complex64_t*, tau,
        PLASMA_Complex64_t*, work,
        double*,             rwork,
        PLASMA_sequence*,    sequence,
        PLASMA_request*,     request);

    return PLASMA_SUCCESS;
}
