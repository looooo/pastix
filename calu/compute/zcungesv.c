/**
 *
 * @file zcungesv.c
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
#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define PLASMA_zlag2c(_descA, _descSB)		      \
  plasma_parallel_call_4(plasma_pzlag2c,	      \
                         PLASMA_desc,      (_descA),  \
                         PLASMA_desc,      (_descSB), \
                         PLASMA_sequence*, sequence,  \
                         PLASMA_request*,  request)

#define PLASMA_clag2z(_descSA, _descB)		      \
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
 *  PLASMA_zcungesv - Solves overdetermined or underdetermined linear systems involving an M-by-N
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
 *  PLASMA_zcungesv first attempts to factorize the matrix in COMPLEX and use this
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
 *  We follow Bjorck's algorithm proposed in "Iterative Refinement of Linear
 *  Least Squares solutions I", BIT, 7:257-278, 1967.4
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Intended usage:
 *          = PlasmaNoTrans:   the linear system involves A;
 *          = PlasmaConjTrans: the linear system involves A**H.
 *          Currently only PlasmaNoTrans is supported.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrices B and X.
 *          NRHS >= 0.
 *
 * @param[in] A
  *          The M-by-N matrix A. This matrix is not modified.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          The M-by-NRHS matrix B of right hand side vectors, stored columnwise. Not modified.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= MAX(1,M,N).
 *
 * @param[out] X
 *          If return value = 0, the solution vectors, stored columnwise.
 *          if M >= N, rows 1 to N of B contain the least squares solution vectors; the residual
 *          sum of squares for the solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column;
 *          if M < N, rows 1 to N of B contain the minimum norm solution vectors;
 *
 * @param[in] LDX
 *          The leading dimension of the array B. LDB >= MAX(1,M,N).
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zcungesv_Tile
 * @sa PLASMA_zcungesv_Tile_Async
 * @sa PLASMA_dsungesv
 * @sa PLASMA_zgels
 *
 ******************************************************************************/
int PLASMA_zcungesv(PLASMA_enum trans, int N, int NRHS,
                    PLASMA_Complex64_t *A, int LDA,
                    PLASMA_Complex64_t *B, int LDB,
                    PLASMA_Complex64_t *X, int LDX, int *ITER)
{
    int NB;
    int status;
    PLASMA_desc  descA;
    PLASMA_desc  descB;
    PLASMA_desc *descT;
    PLASMA_desc  descX;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcungesv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    *ITER = 0;

    /* Check input arguments */
    if (trans != PlasmaNoTrans   && 
        trans != PlasmaConjTrans &&
        trans != PlasmaTrans ) 
    {
        plasma_error("PLASMA_zcungesv", "illegal value of trans");
        return -1;
    }
    if (trans != PlasmaNoTrans) {
        plasma_error("PLASMA_zcungesv", "only PlasmaNoTrans supported");
        return PLASMA_ERR_NOT_SUPPORTED;
    }
    if (N < 0) {
        plasma_error("PLASMA_zcungesv", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        plasma_error("PLASMA_zcungesv", "illegal value of NRHS");
        return -3;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zcungesv", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, N)) {
        plasma_error("PLASMA_zcungesv", "illegal value of LDB");
        return -7;
    }
    if (LDX < max(1, N)) {
        plasma_error("PLASMA_zcungesv", "illegal value of LDX");
        return -9;
    }

    /* Quick return */
    if ( N == 0 )
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZCGELS, N, N, NRHS);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcungesv", "plasma_tune() failed");
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

    /* Allocate workspace */
    PLASMA_Alloc_Workspace_zgels_Tile(N, N, &descT);

    /* Call the native interface */
    status = PLASMA_zcungesv_Tile_Async(PlasmaNoTrans, &descA, descT, &descB, &descX, ITER,
                                        sequence, &request);

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

    PLASMA_Dealloc_Handle_Tile(&descT);
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zcungesv_Tile - Solves symmetric linear system of equations using the tile QR
 *  or the tile LQ factorization and mixed-precision iterative refinement.
 *  Tile equivalent of PLASMA_zcungesv().
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
 *          - If the iterative refinement converged, A is not modified;
 *          - otherwise, it fell back to double precision solution, and
 *          on exit the M-by-N matrix A contains:
 *          if M >= N, A is overwritten by details of its QR factorization as returned by
 *                     PLASMA_zgeqrf;
 *          if M < N, A is overwritten by details of its LQ factorization as returned by
 *                      PLASMA_zgelqf.
 *
 * @param[out] T
 *          On exit:
 *          - if the iterative refinement converged, T is not modified;
 *          - otherwise, it fell back to double precision solution,
 *          and then T is an auxiliary factorization data.
 *
 * @param[in,out] B
 *          On entry, the M-by-NRHS matrix B of right hand side vectors, stored columnwise;
 * @param[in] B
 *          The N-by-NRHS matrix of right hand side matrix B.
 *
 * @param[out] X
 *          If return value = 0, X is the solution vectors, stored columnwise:
 *          if M >= N, rows 1 to N of X contain the least squares solution vectors; the residual
 *          sum of squares for the solution in each column is given by the sum of squares of the
 *          modulus of elements N+1 to M in that column;
 *          if M < N, rows 1 to N of X contain the minimum norm solution vectors;
 *
 * @param[out] ITER
 *          The number of the current iteration in the iterative refinement process
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zcungesv
 * @sa PLASMA_zcungesv_Tile_Async
 * @sa PLASMA_dsungesv_Tile
 * @sa PLASMA_zgels_Tile
 *
 ******************************************************************************/
int PLASMA_zcungesv_Tile(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T,
                         PLASMA_desc *B, PLASMA_desc *X, int *ITER)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcungesv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    status = PLASMA_zcungesv_Tile_Async(trans, A, T, B, X, ITER, sequence, &request);
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
 *  PLASMA_zcungesv_Tile_Async - Solves symmetric linear system of equations using
 *  the tile QR or the tile LQ factorization and mixed-precision iterative refinement.
 *  Non-blocking equivalent of PLASMA_zcungesv_Tile().
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
 * @sa PLASMA_zcungesv
 * @sa PLASMA_zcungesv_Tile
 * @sa PLASMA_dsungesv_Tile_Async
 * @sa PLASMA_zgels_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zcungesv_Tile_Async(PLASMA_enum trans, PLASMA_desc *A, PLASMA_desc *T,
                               PLASMA_desc *B, PLASMA_desc *X, int *ITER,
                               PLASMA_sequence *sequence, PLASMA_request *request)
{
    int N, NB, IB;
    PLASMA_desc descA;
    PLASMA_desc descT;
    PLASMA_desc descB;
    PLASMA_desc descX;
    PLASMA_desc descR, descSA, descST, descSX;
    plasma_context_t *plasma;
    double *work;

    const int    itermax = 30;
    const double bwdmax  = 1.0;
    const PLASMA_Complex64_t negone = -1.0;
    const PLASMA_Complex64_t one = 1.0;
    int iiter;
    double Anorm, cte, eps, Rnorm, Xnorm;
    *ITER=0;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zcungesv_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zcungesv_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zcungesv_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);
    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcungesv_Tile", "invalid first descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descA = *A;
    }
    if (plasma_desc_check(T) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcungesv_Tile", "invalid second descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descT = *T;
    }
    if (plasma_desc_check(B) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcungesv_Tile", "invalid third descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descB = *B;
    }
    if (plasma_desc_check(X) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zcungesv_Tile", "invalid fourth descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    } else {
        descX = *X;
    }
    /* Check input arguments */
    if ( (descA.nb != descA.mb) || (descB.nb != descB.mb) || (descX.nb != descX.mb) ||
         (descA.mb != descB.mb) || (descB.mb != descX.mb) ) {
        plasma_error("PLASMA_zcungesv_Tile", "only square tiles of same size are supported");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    if (trans != PlasmaNoTrans) {
        plasma_error("PLASMA_zcungesv_Tile", "only PlasmaNoTrans supported");
        return PLASMA_ERR_NOT_SUPPORTED;
    }

    /* Set N, NRHS, NB */
    N  = descA.m;
    NB = descA.nb;
    IB = descT.mb;

    work = (double *)plasma_shared_alloc(plasma, PLASMA_SIZE, PlasmaRealDouble);
    if (work == NULL) {
        plasma_error("PLASMA_zcungesv", "plasma_shared_alloc() failed");
        plasma_shared_free(plasma, work);
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    plasma_zdesc_alloc( descR,  NB, NB, descB.m, descB.n, 0, 0, descB.m, descB.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR) );
    plasma_cdesc_alloc( descSA, NB, NB, descA.m, descA.n, 0, 0, descA.m, descA.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR); plasma_desc_mat_free(&descSA) );
    plasma_cdesc_alloc( descST, IB, NB, descT.m, descT.n, 0, 0, descT.m, descT.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR); plasma_desc_mat_free(&descSA); plasma_desc_mat_free(&descST) );
    plasma_cdesc_alloc( descSX, NB, NB, descX.m, descX.n, 0, 0, descX.m, descX.n, plasma_shared_free( plasma, work ); plasma_desc_mat_free(&descR); plasma_desc_mat_free(&descSA); plasma_desc_mat_free(&descST); plasma_desc_mat_free(&descSX) );

    /* Compute some constants */
    PLASMA_zlange(PlasmaInfNorm, descA, Anorm, work);
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

    /* Compute the QR factorization of SA */
    plasma_parallel_call_4(plasma_pcgeqrf,
        PLASMA_desc, descSA,
        PLASMA_desc, descST,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Compute the solve in simple */
    plasma_parallel_call_7(plasma_pcunmqr,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaConjTrans,
        PLASMA_desc, descSA,
        PLASMA_desc, descSX,
        PLASMA_desc, descST,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pctrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex32_t, 1.0,
        PLASMA_desc, descSA,
        PLASMA_desc, descSX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* Convert SX back to double precision */
    PLASMA_clag2z(descSX, descX);

    /* Compute R = B - AX. */
    PLASMA_zlacpy(descB, descR);

    plasma_parallel_call_9(plasma_pzgemm,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNoTrans,
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
        plasma_desc_mat_free(&descST);
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

        plasma_parallel_call_7(plasma_pcunmqr,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaConjTrans,
            PLASMA_desc, descSA,
            PLASMA_desc, descSX,
            PLASMA_desc, descST,
            PLASMA_sequence*, sequence,
            PLASMA_request*, request);

        plasma_parallel_call_9(plasma_pctrsm,
            PLASMA_enum, PlasmaLeft,
            PLASMA_enum, PlasmaUpper,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNonUnit,
            PLASMA_Complex32_t, (PLASMA_Complex32_t)1.0,
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
        plasma_parallel_call_9(plasma_pzgemm,
            PLASMA_enum, PlasmaNoTrans,
            PLASMA_enum, PlasmaNoTrans,
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
            plasma_desc_mat_free(&descST);
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
    plasma_desc_mat_free(&descST);
    plasma_desc_mat_free(&descSX);
    plasma_desc_mat_free(&descR);
    plasma_shared_free(plasma, work);

    /* Single-precision iterative refinement failed to converge to a
       satisfactory solution, so we restart to double precision. */
    PLASMA_zlacpy(descB, descX);

    plasma_parallel_call_4(plasma_pzgeqrf,
        PLASMA_desc, descA,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_7(plasma_pzunmqr,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaConjTrans,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        PLASMA_desc, descT,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_parallel_call_9(plasma_pztrsm,
        PLASMA_enum, PlasmaLeft,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_Complex64_t, (PLASMA_Complex64_t)1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descX,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    return PLASMA_SUCCESS;
}
