/**
 *
 * @file zcposv.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Emmanuel Agullo
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define PLASMA_zlag2c(_descA, _descSB)		      \
  plasma_parallel_call_4(plasma_pzlag2c,	      \
                         PLASMA_desc,      (_descA),  \
                         PLASMA_desc,      (_descSB), \
                         PLASMA_sequence*, sequence,  \
                         PLASMA_request*,  request)

#define PLASMA_clag2z(_descSA, _descB)                \
  plasma_parallel_call_4(plasma_pclag2z,              \
                         PLASMA_desc,      (_descSA), \
                         PLASMA_desc,      (_descB),  \
                         PLASMA_sequence*, sequence,  \
                         PLASMA_request*,  request)

#define PLASMA_zlange(_norm, _descA, _result, _work)   \
  _result = 0;                                         \
  plasma_parallel_call_6(plasma_pzlange,               \
                         PLASMA_enum,      (_norm),    \
                         PLASMA_desc,      (_descA),   \
                         double*,          (_work),    \
                         double*,          &(_result), \
                         PLASMA_sequence*, sequence,   \
                         PLASMA_request*,  request);

#define PLASMA_zlanhe(_norm, _uplo, _descA, _result, _work) \
  _result = 0;                                              \
  plasma_parallel_call_7(plasma_pzlanhe,                    \
                         PLASMA_enum,      (_norm),         \
                         PLASMA_enum,      (_uplo),         \
                         PLASMA_desc,      (_descA),        \
                         double*,          (_work),         \
                         double*,          &(_result),      \
                         PLASMA_sequence*, sequence,        \
                         PLASMA_request*,  request);

#define PLASMA_zlacpy(_descA, _descB)                        \
  plasma_parallel_call_5(plasma_pzlacpy,                     \
                         PLASMA_enum,      PlasmaUpperLower, \
                         PLASMA_desc,      (_descA),         \
                         PLASMA_desc,      (_descB),         \
                         PLASMA_sequence*, sequence,         \
                         PLASMA_request*,  request)

#define PLASMA_zgeadd(_alpha, _descA, _descB)        	\
  plasma_parallel_call_5(plasma_pzgeadd,                	\
                         PLASMA_Complex64_t, (_alpha),	\
                         PLASMA_desc,        (_descA),	\
                         PLASMA_desc,        (_descB),	\
                         PLASMA_sequence*,   sequence,	\
                         PLASMA_request*,    request)

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zcposv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N symmetric positive definite (or Hermitian positive definite
 *  in the complex case) matrix and X and B are N-by-NRHS matrices.
 *  The Cholesky decomposition is used to factor A as
 *
 *    A = U**H * U, if uplo = PlasmaUpper, or
 *    A = L * L**H, if uplo =  PlasmaLower,
 *
 *  where U is an upper triangular matrix and  L is a lower triangular matrix.
 *  The factored form of A is then used to solve the system of equations A * X = B.
 *
 *  PLASMA_zcposv first attempts to factorize the matrix in COMPLEX and use this
 *  factorization within an iterative refinement procedure to produce a
 *  solution with COMPLEX*16 normwise backward error quality (see below).
 *  If the approach fails the method switches to a COMPLEX*16
 *  factorization and solve.
 *
 *  The iterative refinement is not going to be a winning strategy if
 *  the ratio COMPLEX performance over COMPLEX*16 performance is too
 *  small. A reasonable strategy should take the number of right-hand
 *  sides and the size of the matrix into account. This might be done
 *  with a call to ILAENV in the future. Up to now, we always try
 *  iterative refinement.
 *
 *  The iterative refinement process is stopped if ITER > ITERMAX or
 *  for all the RHS we have: RNRM < N*XNRM*ANRM*EPS*BWDMAX
 *  where:
 *
 *  - ITER is the number of the current iteration in the iterative refinement process
 *  - RNRM is the infinity-norm of the residual
 *  - XNRM is the infinity-norm of the solution
 *  - ANRM is the infinity-operator-norm of the matrix A
 *  - EPS is the machine epsilon returned by DLAMCH('Epsilon').
 *
 *  Actually, in its current state (PLASMA 2.1.0), the test is slightly relaxed.
 *
 *  The values ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 respectively.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in] A
 *          The N-by-N symmetric positive definite (or Hermitian) coefficient matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          This matrix is not modified.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] B
 *          The N-by-NRHS matrix of right hand side matrix B.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[out] X
 *          If return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDX
 *          The leading dimension of the array B. LDX >= max(1,N).
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zcposv_Tile
 * @sa PLASMA_zcposv_Tile_Async
 * @sa PLASMA_dsposv
 * @sa PLASMA_zposv
 *
 ******************************************************************************/
int PLASMA_zcposv(PLASMA_enum uplo, int N, int NRHS,
                  PLASMA_Complex64_t *A, int LDA,
                  PLASMA_Complex64_t *B, int LDB,
                  PLASMA_Complex64_t *X, int LDX, int *ITER)
{
    int NB;
    int status;
    PLASMA_desc  descA;
    PLASMA_desc  descB;
    PLASMA_desc  descX;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcposv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    *ITER = 0;

    /* Check input arguments */
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_zcposv", "illegal value of uplo");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zcposv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_zcposv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zcposv", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_zcposv", "illegal value of LDB");
        return -7;
    }
    if (LDX < max(1, N)) {
        plasma_error("PLASMA_zcposv", "illegal value of LDX");
        return -9;
    }
    /* Quick return - currently NOT equivalent to LAPACK's
     * LAPACK does not have such check for ZCPOSV */
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZCPOSV, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcposv", "plasma_tune() failed");
        return status;
    }

    NB = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    /* DOUBLE PRECISION INITIALIZATION */
    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,    sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
        plasma_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                             plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)) );
        plasma_zdesc_alloc(  descX, NB, NB, N, NRHS, 0, 0, N, NRHS, plasma_desc_mat_free(&(descA)); plasma_desc_mat_free(&(descB)); plasma_desc_mat_free(&(descX)) );
    } else {
        plasma_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,   
                            sequence, &request);
        plasma_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,
                            sequence, &request);

        descX = plasma_desc_init(
            PlasmaComplexDouble, NB, NB, (NB*NB), 
            LDX, NRHS, 0, 0, N, NRHS);
        descX.mat = X;
    }

    /* Call the native interface */
    status = PLASMA_zcposv_Tile_Async(uplo, &descA, &descB, &descX, ITER, sequence, &request);

    if (status == PLASMA_SUCCESS) {
        if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
            plasma_zooptile2lap( descX, X, NB, NB, LDX, NRHS,  sequence, &request);
            plasma_dynamic_sync();
            plasma_desc_mat_free(&descA);
            plasma_desc_mat_free(&descB);
            plasma_desc_mat_free(&descX);
        } else {
            plasma_ziptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);
            plasma_ziptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);
            plasma_ziptile2lap( descX, X, NB, NB, LDX, NRHS,  sequence, &request);
            plasma_dynamic_sync();
        }
    }

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}


/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zcposv_Tile - Solves a symmetric positive definite or Hermitian positive definite
 *  system of linear equations using the Cholesky factorization and mixed-precision iterative refinement.
 *  Tile equivalent of PLASMA_zcposv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = PlasmaUpper: Upper triangle of A is stored;
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in,out] A
 *          On entry, the N-by-N symmetric positive definite (or Hermitian) coefficient matrix A.
 *          If uplo = PlasmaUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          - If the iterative refinement converged, A is not modified;
 *          - otherwise, it falled backed to double precision solution,
 *
 * @param[in] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *
 * @param[out] X
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zcposv
 * @sa PLASMA_zcposv_Tile_Async
 * @sa PLASMA_dsposv_Tile
 * @sa PLASMA_zposv_Tile
 *
 ******************************************************************************/
int PLASMA_zcposv_Tile(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B,
                       PLASMA_desc *X, int *ITER)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcposv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    status = PLASMA_zcposv_Tile_Async(uplo, A, B, X, ITER, sequence, &request);
    if (status != PLASMA_SUCCESS)
        return status;
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zcposv_Tile_Async - Solves a symmetric positive definite or Hermitian
 *  positive definite system of linear equations using the Cholesky factorization
 *  and mixed-precision iterative refinement.
 *  Non-blocking equivalent of PLASMA_zcposv_Tile().
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
 * @sa PLASMA_zcposv
 * @sa PLASMA_zcposv_Tile
 * @sa PLASMA_dsposv_Tile_Async
 * @sa PLASMA_zposv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zcposv_Tile_Async(PLASMA_enum uplo, PLASMA_desc *A, PLASMA_desc *B,
                             PLASMA_desc *X, int *ITER,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    int N, NB;
    PLASMA_desc descA;
    PLASMA_desc descB;
    PLASMA_desc descX;
    plasma_context_t *plasma;
    double *work;
    PLASMA_desc descR, descSA, descSX;

    const int itermax = 30;
    const double bwdmax = 1.0;
    const PLASMA_Complex64_t negone = -1.0;
    const PLASMA_Complex64_t one = 1.0;
    int iiter;
    double Anorm, cte, eps, Rnorm, Xnorm;
    *ITER=0;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcposv_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zcposv_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zcposv_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcposv_Tile_Async", "invalid first descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descA = *A;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcposv_Tile_Async", "invalid second descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descB = *B;
    }
    if (plasma_desc_check(X) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcposv_Tile_Async", "invalid third descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descX = *X;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb || descB.nb != descB.mb || descX.nb != descX.mb) {
        plasma_error("PLASMA_zcposv_Tile_Async", "only square tiles supported");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("PLASMA_zcposv_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return - currently NOT equivalent to LAPACK's
     * LAPACK does not have such check for DPOSV */

/*
    if (min(N, NRHS) == 0)
        return PLASMA_SUCCESS;
*/

    /* Set N, NRHS */
    N  = descA.m;
    NB = descA.nb;

    work = (double *)plasma_shared_alloc(plasma, PLASMA_SIZE, PlasmaRealDouble);
    if (work == NULL) {
        plasma_error("PLASMA_zcposv_Tile_Async", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, work);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    plasma_zdesc_alloc( descR,  NB, NB, descB.m, descB.n, 0, 0, descB.m, descB.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR) );
    plasma_cdesc_alloc( descSA, NB, NB, descA.m, descA.n, 0, 0, descA.m, descA.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR); plasma_desc_mat_free(&descSA) );
    plasma_cdesc_alloc( descSX, NB, NB, descX.m, descX.n, 0, 0, descX.m, descX.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR); plasma_desc_mat_free(&descSA); plasma_desc_mat_free(&descSX) );

    /* Compute some constants */
    PLASMA_zlanhe(PlasmaInfNorm, uplo, descA, Anorm, work);
    eps = LAPACKE_dlamch_work('e');

    /* Convert B from double precision to single precision and store
       the result in SX. */
    PLASMA_zlag2c(descB, descSX);
    if (sequence->status != PLASMA_SUCCESS)
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Convert A from double precision to single precision and store
       the result in SA. */
    PLASMA_zlag2c(descA, descSA);
    if (sequence->status != PLASMA_SUCCESS)
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Compute the Cholesky factorization of SA */
    plasma_parallel_call_4(plasma_pcpotrf,
        PLASMA_enum, uplo,
        PLASMA_desc, descSA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Solve the system SA*SX = SB */
    /* Forward substitution */
    plasma_parallel_call_9(plasma_pctrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaConjTrans : PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex32_t, 1.0,
        PLASMA_desc, descSA,
        PLASMA_desc, descSX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Backward substitution */
    plasma_parallel_call_9(plasma_pctrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaConjTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex32_t, 1.0,
        PLASMA_desc, descSA,
        PLASMA_desc, descSX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Convert SX back to double precision */
    PLASMA_clag2z(descSX, descX);

    /* Compute R = B - AX. */
    PLASMA_zlacpy(descB,descR);
    plasma_parallel_call_9(plasma_pzhemm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_Complex64_t, negone,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        PLASMA_Complex64_t, one,
        PLASMA_desc, descR,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Check whether the NRHS normwise backward error satisfies the
       stopping criterion. If yes return. Note that ITER=0 (already set). */
    PLASMA_zlange(PlasmaInfNorm, descX, Xnorm, work);
    PLASMA_zlange(PlasmaInfNorm, descR, Rnorm, work);

    /* Wait the end of Anorm, Xnorm and Bnorm computations */
    plasma_dynamic_sync();

    cte = Anorm*eps*((double) N)*bwdmax;
    if (Rnorm < Xnorm * cte){
        /* The NRHS normwise backward errors satisfy the
           stopping criterion. We are good to exit. */
        plasma_desc_mat_free(&descSA);
        plasma_desc_mat_free(&descSX);
        plasma_desc_mat_free(&descR);
        plasma_shared_free(plasma, work);
        return PLASMA_SUCCESS;
    }

    /* Iterative refinement */
    for (iiter = 0; iiter < itermax; iiter++){

        /* Convert R from double precision to single precision
           and store the result in SX. */
        PLASMA_zlag2c(descR, descSX);

        /* Solve the system SA*SX = SR */
        /* Forward substitution */
        plasma_parallel_call_9(plasma_pctrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, uplo,
            PLASMA_enum, uplo == PlasmaUpper ? PlasmaConjTrans : PlasmaNoTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex32_t, 1.0,
            PLASMA_desc, descSA,
            PLASMA_desc, descSX,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Backward substitution */
        plasma_parallel_call_9(plasma_pctrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, uplo,
            PLASMA_enum, uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaConjTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex32_t, 1.0,
            PLASMA_desc, descSA,
            PLASMA_desc, descSX,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Convert SX back to double precision and update the current
           iterate. */
        PLASMA_clag2z(descSX, descR);
        PLASMA_zgeadd(one, descR, descX);

        /* Compute R = B - AX. */
        PLASMA_zlacpy(descB,descR);
        plasma_parallel_call_9(plasma_pzhemm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, uplo,
            PLASMA_Complex64_t, negone,
            PLASMA_desc, descA,
            PLASMA_desc, descX,
            PLASMA_Complex64_t, one,
            PLASMA_desc, descR,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        /* Check whether the NRHS normwise backward errors satisfy the
           stopping criterion. If yes, set ITER=IITER>0 and return. */
        PLASMA_zlange(PlasmaInfNorm, descX, Xnorm, work);
        PLASMA_zlange(PlasmaInfNorm, descR, Rnorm, work);

        /* Wait the end of Xnorm and Bnorm computations */
        plasma_dynamic_sync();

        if (Rnorm < Xnorm * cte){
            /* The NRHS normwise backward errors satisfy the
               stopping criterion. We are good to exit. */
            *ITER = iiter;

            plasma_desc_mat_free(&descSA);
            plasma_desc_mat_free(&descSX);
            plasma_desc_mat_free(&descR);
            plasma_shared_free(plasma, work);
            return PLASMA_SUCCESS;
        }
    }

    /* We have performed ITER=itermax iterations and never satisified
       the stopping criterion, set up the ITER flag accordingly and
       follow up on double precision routine. */
    *ITER = -itermax - 1;

    plasma_desc_mat_free(&descSA);
    plasma_desc_mat_free(&descSX);
    plasma_desc_mat_free(&descR);
    plasma_shared_free(plasma, work);

    /* Single-precision iterative refinement failed to converge to a
       satisfactory solution, so we resort to double precision. */

    plasma_parallel_call_4(plasma_pzpotrf,
        PLASMA_enum, uplo,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    PLASMA_zlacpy(descB,descX);

    plasma_parallel_call_9(plasma_pztrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaConjTrans : PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex64_t, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pztrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, uplo,
        PLASMA_enum, uplo == PlasmaUpper ? PlasmaNoTrans : PlasmaConjTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex64_t, 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
