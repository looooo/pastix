/**
 * @file zhegvd.c
 *
 *  PLASMA computational routines
 *  Release Date: November, 15th 2009
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
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
 *  PLASMA_zhegvd - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also positive
 *  definite.
 *
 *******************************************************************************
 *
 * @param[in] itype
 *          Intended usage:
 *          = 1: A*x=(lambda)*B*x
 *          = 2: A*Bx=(lambda)*x
 *          = 3: B*A*x=(lambda)*x
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A and B are stored;
 *          = PlasmaLower: Lower triangle of A and B are stored.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
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
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the matrix Z of eigenvectors.
 *          The eigenvectors are normalized as follows:
 *          if ITYPE = 1 or 2, Z**H*B*Z = I;
 *          if ITYPE = 3,      Z**H*inv(B)*Z = I.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in,out] B
 *          On entry, the symmetric (or Hermitian) positive definite
 *          matrix B.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of B contains the upper triangular part of the matrix
 *          B, and the strictly lower triangular part of B is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of B contains the lower triangular part of the matrix
 *          B, and the strictly upper triangular part of B is not
 *          referenced.
 *          On exit, if return value <= N, the part of B containing
 *          the matrix is overwritten by the triangular factor U or L
 *          from the Cholesky factorization B = U**H*U or B = L*L**H.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDA >= max(1,N).
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[in, out] descT
 *          On entry, descriptor as return by PLASMA_Alloc_Workspace_zhegvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 * @param[in] LDQ
 *          The leading dimension of Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval <=N if INFO = i, plasma_zhegvd failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *          \retval >N if INFO = N + i, for 1 <= i <= N, then the leading
 *                     minor of order i of B is not positive definite.
 *                     The factorization of B could not be completed and
 *                     no eigenvalues or eigenvectors were computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zhegvd_Tile
 * @sa PLASMA_zhegvd_Tile_Async
 * @sa PLASMA_chegvd
 * @sa PLASMA_dsygvd
 * @sa PLASMA_ssygvd
 *
 ******************************************************************************/
int PLASMA_zhegvd(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB,
                 double *W,
                 PLASMA_desc *descT,
                 PLASMA_Complex64_t *Q, int LDQ)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA, descB, descQ;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA_zhegvd", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (itype != 1 && itype != 2 && itype != 3) {
        plasma_error("PLASMA_zhegvd", "Illegal value of itype");
        return -1;
    }
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zhegvd", "illegal value of jobz");
        return -2;
    }
    if (uplo != PlasmaLower && uplo!= PlasmaUpper) {
        plasma_error("PLASMA_zhegvd", "only PlasmaLower supported");
        return -3;
    }
    if (N < 0) {
        plasma_error("PLASMA_zhegvd", "illegal value of N");
        return -4;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zhegvd", "illegal value of LDA");
        return -6;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_zhegvd", "illegal value of LDB");
        return -8;
    }
    if (LDQ < max(1, N)) {
        plasma_error("PLASMA_zhegvd", "illegal value of LDQ");
        return -12;
    }

    /* Quick return */
    if (N == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on N; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZHEGV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegvd", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)) );
        /* No conversion for Q */
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Q is used only as an output */
    descQ = plasma_desc_init(
        PlasmaComplexDouble, NB, NB, NB*NB,
        LDQ, N, 0, 0, N, N);
    descQ.mat = Q;

    /* Call the tile interface */
    PLASMA_zhegvd_Tile_Async(itype, jobz, uplo,
                            &descA, &descB, W,
                            descT, &descQ,
                            sequence, &request);

    plasma_ziptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);
    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N, sequence, &request);
        plasma_zooptile2lap( descB, B, NB, NB, LDB, N, sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
        plasma_desc_mat_free(&descB);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N, sequence, &request);
        plasma_ziptile2lap( descB, B, NB, NB, LDB, N, sequence, &request);
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
 *  PLASMA_zhegvd_Tile - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also
 *  positive definite.
 *
 *  Tile equivalent of PLASMA_zhegvd().
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
 *
 * @param[in] jobz
 *          Intended usage:
 *          = PlasmaNoVec: computes eigenvalues only;
 *          = PlasmaVec: computes eigenvalues and eigenvectors.
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or
 *          lower triangular:
 *          = PlasmaUpper: Upper triangle of A and B are stored;
 *          = PlasmaLower: Lower triangle of A and B are stored.
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
 *          On exit, if jobz = PlasmaVec, then if return value = 0, A
 *          contains the matrix Z of eigenvectors.
 *          The eigenvectors are normalized as follows:
 *          if ITYPE = 1 or 2, Z**H*B*Z = I;
 *          if ITYPE = 3,      Z**H*inv(B)*Z = I.
 *          If jobz = PlasmaNoVec, then on exit the lower triangle (if
 *          uplo = PlasmaLower) or the upper triangle (if uplo =
 *          PlasmaUpper) of A, including the diagonal, is destroyed.
 *
 * @param[in,out] B
 *          On entry, the symmetric (or Hermitian) positive definite
 *          matrix B.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular
 *          part of B contains the upper triangular part of the matrix
 *          B, and the strictly lower triangular part of B is not
 *          referenced.
 *          If uplo = PlasmaLower, the leading N-by-N lower triangular
 *          part of B contains the lower triangular part of the matrix
 *          B, and the strictly upper triangular part of B is not
 *          referenced.
 *          On exit, if return value <= N, the part of B containing
 *          the matrix is overwritten by the triangular factor U or L
 *          from the Cholesky factorization B = U**H*U or B = L*L**H.
 *
 * @param[in,out] T
 *          On entry, descriptor as return by
 *          PLASMA_Alloc_Workspace_zhegvd
 *          On exit, contains auxiliary factorization data.
 *
 * @param[out] W
 *          On exit, if info = 0, the eigenvalues.
 *
 * @param[out] Q
 *          On exit, if jobz = PlasmaVec and info = 0, the eigenvectors.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval <=N if INFO = i, plasma_zhegvd failed to converge; i
 *               off-diagonal elements of an intermediate tridiagonal
 *               form did not converge to zero.
 *          \retval >N if INFO = N + i, for 1 <= i <= N, then the leading
 *                     minor of order i of B is not positive definite.
 *                     The factorization of B could not be completed and
 *                     no eigenvalues or eigenvectors were computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zhegvd
 * @sa PLASMA_zhegvd_Tile_Async
 * @sa PLASMA_chegvd_Tile
 * @sa PLASMA_dsygvd_Tile
 * @sa PLASMA_ssygvd_Tile
 *
 ******************************************************************************/
int PLASMA_zhegvd_Tile(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo,
                      PLASMA_desc *A, PLASMA_desc *B, double *W,
                      PLASMA_desc *T, PLASMA_desc *Q)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhegvd_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zhegvd_Tile_Async(itype, jobz, uplo, A, B, W, T, Q, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zhegvd_Tile - Computes all eigenvalues and, optionally,
 *  eigenvectors of a complex generalized Hermitian-definite
 *  eigenproblem of the form: A*x=(lambda)*B*x, A*Bx=(lambda)*x, or
 *  B*A*x=(lambda)*x.
 *  Here A and B are assumed to be Hermitian and B is also
 *  positive definite.
 *
 *  Non-blocking equivalent of PLASMA_zhegvd_Tile().
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
 * @sa PLASMA_zhegvd
 * @sa PLASMA_zhegvd_Tile
 * @sa PLASMA_chegvd_Tile_Async
 * @sa PLASMA_dsygvd_Tile_Async
 * @sa PLASMA_ssygvd_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zhegvd_Tile_Async(PLASMA_enum itype, PLASMA_enum jobz, PLASMA_enum uplo,
                             PLASMA_desc *A, PLASMA_desc *B, double *W,
                             PLASMA_desc *T, PLASMA_desc *Q,
PLASMA_sequence *sequence, PLASMA_request *request)
{
    int status;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descQ;
    PLASMA_enum trans;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zhegvd_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zhegvd_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zhegvd_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "invalid second descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descB = *B;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "invalid T descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (plasma_desc_check(Q) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "invalid Q descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descQ = *Q;
    }

    /* Check input arguments */
    if (itype != 1 && itype != 2 && itype != 3) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "Illegal value of itype");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (jobz != PlasmaNoVec && jobz != PlasmaVec) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "illegal value of jobz");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != PlasmaLower && uplo != PlasmaUpper) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "illegal value of uplo");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zhegvd_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    if ( (jobz == PlasmaVec) && (descA.nb != descA.mb) ) {
        plasma_error("PLASMA_zhegv_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /*
     * Form a Cholesky factorization of B.
     */
    plasma_parallel_call_4(plasma_pzpotrf,
        PLASMA_enum,      uplo,
        PLASMA_desc,      descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    status = sequence->status;
    if (status != 0){
       status = descA.n + status;
       return status;
    }

    /*
     * Transform problem to standard eigenvalue problem and solve.
     */
    plasma_dynamic_call_6(plasma_pzhegst,
        PLASMA_enum,      itype,
        PLASMA_enum,      uplo,
        PLASMA_desc,      descA,
        PLASMA_desc,      descB,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    /*
     * Call Eigensolver
     */
    status = PLASMA_zheevd_Tile_Async( jobz, uplo,
                                       A, W, T,
                                       (PLASMA_Complex64_t*)(Q->mat), Q->lm,
                                       sequence, request );

    /*
     * Backtransform eigenvectors to the original problem.
     */
    if( jobz == PlasmaVec ) {
        if ( status != PLASMA_SUCCESS ) {
            plasma_error("PLASMA_zheegvd", "zhegvd");
            return status;
        }

        /*
         * Convert back the lapack layout used in zheev to tile layout
         */
        plasma_ziplap2tile( descQ, (PLASMA_Complex64_t*)Q->mat, descQ.mb, descQ.nb, descQ.lm, descQ.ln,
                            0, 0, descQ.m, descQ.n, sequence, request );

        if( (itype == 1) || (itype == 2) ) {
            /*
             * For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
             * backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
             */
             if( uplo == PlasmaUpper )
                 trans = PlasmaNoTrans;
             else
                 trans = PlasmaConjTrans;

             plasma_parallel_call_9(plasma_pztrsm,
                 PLASMA_enum,        PlasmaLeft,
                 PLASMA_enum,        uplo,
                 PLASMA_enum,        trans,
                 PLASMA_enum,        PlasmaNonUnit,
                 PLASMA_Complex64_t, 1.0,
                 PLASMA_desc,        descB,
                 PLASMA_desc,        descQ,
                 PLASMA_sequence*,   sequence,
                 PLASMA_request*,    request);
        }
        else if( itype == 3 ) {
            /*
             * For B*A*x=(lambda)*x;
             * backtransform eigenvectors: x = L*y or U**H *y
             */
             if( uplo == PlasmaUpper )
                 trans = PlasmaConjTrans;
             else
                 trans = PlasmaNoTrans;

             plasma_dynamic_call_9(plasma_pztrmm,
                 PLASMA_enum,        PlasmaLeft,
                 PLASMA_enum,        uplo,
                 PLASMA_enum,        trans,
                 PLASMA_enum,        PlasmaNonUnit,
                 PLASMA_Complex64_t, 1.0,
                 PLASMA_desc,        descB,
                 PLASMA_desc,        descQ,
                 PLASMA_sequence*,   sequence,
                 PLASMA_request*,    request);
        }
    }

    return PLASMA_SUCCESS;
}
