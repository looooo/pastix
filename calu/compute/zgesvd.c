/**
 *
 * @file zgesvd.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zgesvd - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *
 *       A = U * SIGMA * transpose(V)
 *
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *
 *  Note that the routine returns V**T, not V.
 *  Not LAPACK Compliant for now!
 *  Note: Only PlasmaNoVec supported!
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *          Note: Only PlasmaNoVec supported!
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**H.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *          Note: Only PlasmaNoVec supported!
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
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**H (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] S
 *          The double precision singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[in] LDU
 *          The leading dimension of the array U.  LDU >= 1; if
 *          JOBU = 'S' or 'A', LDU >= M.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**H;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**H (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[in] LDVT
 *         The leading dimension of the array VT.  LDVT >= 1; if
 *         JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_zgesvd
 *          On exit, contains auxiliary factorization data.
 *
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgesvd_Tile
 * @sa PLASMA_zgesvd_Tile_Async
 * @sa PLASMA_cgesvd
 * @sa PLASMA_dgesvd
 * @sa PLASMA_sgesvd
 *
 ******************************************************************************/
int PLASMA_zgesvd(PLASMA_enum jobu, PLASMA_enum jobvt, int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  double *S,
                  PLASMA_Complex64_t *U, int LDU,
                  PLASMA_Complex64_t *VT, int LDVT,
                  PLASMA_desc *descT)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descU, descVT;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgesvd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (jobu != PlasmaNoVec  && jobu !=PlasmaVec) {
        plasma_error("PLASMA_zgesvd", "illegal value of jobu");
        return -1;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_zgesvd", "illegal value of jobvt");
        return -2;
    }
    if (M < 0) {
        plasma_error("PLASMA_zgesvd", "illegal value of M");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgesvd", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgesvd", "illegal value of LDA");
        return -6;
    }
    if (LDU < 1) {
        plasma_error("PLASMA_zgesvd", "illegal value of LDU");
        return -9;
    }
    if (LDVT < 1) {
        plasma_error("PLASMA_zgesvd", "illegal value of LDVT");
        return -11;
    }
    /* Quick return */
    if (min(M, N) == 0) {
        return PLASMA_SUCCESS;
    }

    if (jobu == PlasmaVec) {
        plasma_error("PLASMA_zgesvd", "computing the singular vectors is not supported in this version");
        return -1;
    }
    if (jobvt == PlasmaVec) {
        plasma_error("PLASMA_zgesvd", "computing the singular vectors is not supported in this version");
        return -2;
    }

    /* Tune NB & IB depending on M & N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGESVD, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgesvd", "plasma_tune() failed");
        return status;
    }

    /* Set MT, NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        if (jobu == PlasmaVec){
            plasma_zooplap2tile( descU,   U, NB, NB,  LDU, M, 0, 0, M, M, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descU)));
        }
        if (jobvt == PlasmaVec){
            plasma_zooplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descU)); plasma_desc_mat_free(&(descVT)));
        }
    } else {
        plasma_ziplap2tile( descA,   A, NB, NB,  LDA, N, 0, 0, M, N,
                            sequence, &request);
        if (jobu == PlasmaVec){
            plasma_ziplap2tile( descU,   U, NB, NB,  LDU, M, 0, 0, M, M,
                            sequence, &request);
        }
        if (jobvt == PlasmaVec){
            plasma_ziplap2tile( descVT, VT, NB, NB, LDVT, N, 0, 0, N, N,
                            sequence, &request);
        }
    }

    /* Call the tile interface */
    PLASMA_zgesvd_Tile_Async(jobu, jobvt, &descA, S, &descU, &descVT, descT, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA,   A, NB, NB,  LDA, N,  sequence, &request);
        if (jobu == PlasmaVec){
            plasma_zooptile2lap( descU,   U, NB, NB,  LDU, M,  sequence, &request);
        }
        if (jobvt == PlasmaVec){
            plasma_zooptile2lap( descVT, VT, NB, NB, LDVT, N,  sequence, &request);
        }
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        if (jobu == PlasmaVec){
            plasma_desc_mat_free(&descU);
        }
        if (jobvt == PlasmaVec){
            plasma_desc_mat_free(&descVT);
        }
    } else {
        plasma_ziptile2lap( descA,   A, NB, NB,  LDA, N,  sequence, &request);
        if (jobu == PlasmaVec){
            plasma_ziptile2lap( descU,   U, NB, NB,  LDU, M,  sequence, &request);
        }
        if (jobvt == PlasmaVec){
            plasma_ziptile2lap( descVT, VT, NB, NB, LDVT, N,  sequence, &request);
        }
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
 *  PLASMA_zgesvd_Tile - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Tile equivalent of PLASMA_zgesvd().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] jobu
 *          Specifies options for computing all or part of the matrix U.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *
 * @param[in] jobvt
 *          Specifies options for computing all or part of the matrix V**H.
 *          Intended usage:
 *          = PlasmaVec: all M columns of U are returned in array U;
 *          = PlasmaNoVec: no columns of U (no left singular vectors) are
 *                     computed.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit,
 *          if JOBU = 'O',  A is overwritten with the first min(m,n)
 *                          columns of U (the left singular vectors,
 *                          stored columnwise);
 *          if JOBVT = 'O', A is overwritten with the first min(m,n)
 *                          rows of V**H (the right singular vectors,
 *                          stored rowwise);
 *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
 *                          are destroyed.
 *
 * @param[out] S
 *          The singular values of A, sorted so that S(i) >= S(i+1).
 *
 * @param[out] U
 *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
 *          If JOBU = 'A', U contains the M-by-M unitary matrix U;
 *          if JOBU = 'S', U contains the first min(m,n) columns of U
 *          (the left singular vectors, stored columnwise);
 *          if JOBU = 'N' or 'O', U is not referenced.
 *
 * @param[out] VT
 *         If JOBVT = 'A', VT contains the N-by-N unitary matrix
 *         V**H;
 *         if JOBVT = 'S', VT contains the first min(m,n) rows of
 *         V**H (the right singular vectors, stored rowwise);
 *         if JOBVT = 'N' or 'O', VT is not referenced.
 *
 * @param[out] T
 *         On exit, auxiliary factorization data.
 *
 *******************************************************************************
 *
 * @return
 *          \return PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgesvd
 * @sa PLASMA_zgesvd_Tile_Async
 * @sa PLASMA_cgesvd_Tile
 * @sa PLASMA_dgesvd_Tile
 * @sa PLASMA_sgesvd_Tile
 *
 ******************************************************************************/
int PLASMA_zgesvd_Tile(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A,
                      double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgesvd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgesvd_Tile_Async(jobu, jobvt, A, S, U, VT, T, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgesvd_Tile_Async - computes the singular value decomposition (SVD) of a complex
 *  M-by-N matrix A, optionally computing the left and/or right singular
 *  vectors.
 *  Non-blocking equivalent of PLASMA_zgesvd_Tile().
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
 * @sa PLASMA_zgesvd
 * @sa PLASMA_zgesvd_Tile
 * @sa PLASMA_cgesvd_Tile_Async
 * @sa PLASMA_dgesvd_Tile_Async
 * @sa PLASMA_sgesvd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgesvd_Tile_Async(PLASMA_enum jobu, PLASMA_enum jobvt, PLASMA_desc *A,
                             double *S, PLASMA_desc *U, PLASMA_desc *VT, PLASMA_desc *T,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA ;
    PLASMA_desc descT ;
    PLASMA_desc descU, descVT;
    double *E;
    int minMN;
    int NCVT = 0;
    int NRU = 0;
    int NCC = 0;

    plasma_context_t *plasma;
    plasma = plasma_context_self();

    if (jobu != PlasmaNoVec  && jobu != PlasmaVec) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "illegal value of jobu");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (jobvt != PlasmaNoVec && jobvt != PlasmaVec) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "illegal value of jobvt");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgesvd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgesvd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgesvd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if ((jobu != PlasmaNoVec) && (plasma_desc_check(U) != PLASMA_SUCCESS)) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((jobvt != PlasmaNoVec) && (plasma_desc_check(VT) != PLASMA_SUCCESS) ) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "invalid third descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "invalid fourth descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descT = *T;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (( (jobu != PlasmaNoVec) && (U->nb != U->mb) )  || ( (jobvt != PlasmaNoVec) && (VT->nb != VT->mb) )) {
        plasma_error("PLASMA_zgesvd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ((jobu == PlasmaVec) || (jobvt == PlasmaVec) ){
        plasma_error("PLASMA_zgesvd_Tile_Async", "computing the singular vectors is not supported in this version");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    if (jobu == PlasmaVec) {
        descU = *U;
    }
    if (jobvt == PlasmaVec) {
        descVT = *VT;
    }

    minMN = min(descA.m, descA.n);
    E = (double *)plasma_shared_alloc(plasma, minMN-1, PlasmaRealDouble);

    /*
     * 1: Reduction to BAND bidiagonal form
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
     * Build the U of the first stage
     */
    if ( jobu == PlasmaVec ) {
        /* Initialize U to Identity */
        plasma_dynamic_call_6(plasma_pzlaset,
            PLASMA_enum, PlasmaUpperLower,
            PLASMA_Complex64_t, 0.0,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descU,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Accumulate the transformations from the first stage */
        /* if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
        /*     plasma_dynamic_call_6(plasma_pzungbr, */
        /*         PLASMA_enum, PlasmaLeft, */
        /*         PLASMA_desc, descA, */
        /*         PLASMA_desc, descU, */
        /*         PLASMA_desc, descT, */
        /*         PLASMA_sequence*, sequence, */
        /*         PLASMA_request*, request); */
        /* } else { */
        /*     plasma_dynamic_call_6(plasma_pzungbrrh, */
        /*         PLASMA_enum, PlasmaLeft, */
        /*         PLASMA_desc, descA, */
        /*         PLASMA_desc, descU, */
        /*         PLASMA_desc, descT, */
        /*         PLASMA_sequence*, sequence, */
        /*         PLASMA_request*, request); */
        /* } */
    }

    /*
     * Build the VT of the first stage
     */
    if ( jobvt == PlasmaVec ) {
        /* Initialize VT to Identity */
        plasma_dynamic_call_6(plasma_pzlaset,
            PLASMA_enum, PlasmaUpperLower,
            PLASMA_Complex64_t, 0.0,
            PLASMA_Complex64_t, 1.0,
            PLASMA_desc, descVT,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Accumulate the transformations from the first stage */
        /* if (plasma->householder == PLASMA_FLAT_HOUSEHOLDER) { */
        /*     plasma_dynamic_call_6(plasma_pzungbr, */
        /*         PLASMA_enum, PlasmaRight, */
        /*         PLASMA_desc, descA, */
        /*         PLASMA_desc, descVT, */
        /*         PLASMA_desc, descT, */
        /*         PLASMA_sequence*, sequence, */
        /*         PLASMA_request*, request); */
        /* } */
        /* else { */
        /*     plasma_dynamic_call_6(plasma_pzungbrrh, */
        /*         PLASMA_enum, PlasmaRight, */
        /*         PLASMA_desc, descA, */
        /*         PLASMA_desc, descVT, */
        /*         PLASMA_desc, descT, */
        /*         PLASMA_sequence*, sequence, */
        /*         PLASMA_request*, request); */
        /* } */
    }

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
     * 2: Reduction from BAND bidiagonal to the final condensed form
     */
    plasma_dynamic_call_7(plasma_pzgebrd_tb2bd,
        PLASMA_enum, descA.m >= descA.n ? PlasmaUpper : PlasmaLower,
        PLASMA_desc, descA,
        double*, S,
        double*, E,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /*
     * Compute the singular values ONLY for now
     */
    plasma_dynamic_sync();
    if (descA.m >= descA.n)
       LAPACKE_zbdsqr(
              LAPACK_COL_MAJOR, lapack_const(PlasmaUpper),
              minMN, NCVT, NRU, NCC,
              S, E,
              NULL, 1, NULL, 1, NULL, 1 );
    else {
       LAPACKE_zbdsqr(
              LAPACK_COL_MAJOR, lapack_const(PlasmaLower),
              minMN, NCVT, NRU, NCC,
              S, E,
              NULL, 1, NULL, 1, NULL, 1 );
    }

    plasma_shared_free(plasma, E);

    return PLASMA_SUCCESS;
}
