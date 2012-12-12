/**
 *
 * @file zgetri.c
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
 * PLASMA_zgetri - Computes the inverse of a matrix using the LU factorization
 * computed by PLASMA_zgetrf.
 * This method inverts U and then computes inv(A) by solving the system
 * inv(A)*L = inv(U) for inv(A).
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the triangular factor L or U from the
 *          factorization A = P*L*U as computed by PLASMA_zgetrf.
 *          On exit, if return value = 0, the inverse of the original
 *          matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[in] IPIV
 *          The pivot indices that define the permutations
 *          as returned by PLASMA_zgetrf.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, the (i,i) element of the factor U is
 *                exactly zero; The matrix is singular 
 *                and its inverse could not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetri_Tile
 * @sa PLASMA_zgetri_Tile_Async
 * @sa PLASMA_cgetri
 * @sa PLASMA_dgetri
 * @sa PLASMA_sgetri
 * @sa PLASMA_zgetrf
 *
 ******************************************************************************/
int PLASMA_zgetri(int N,
                  PLASMA_Complex64_t *A, int LDA,
                  int *IPIV)
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;
    PLASMA_desc descW;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetri", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        plasma_error("PLASMA_zgetri", "illegal value of N");
        return -1;
    }
    if (LDA < max(1, N)) {
        plasma_error("PLASMA_zgetri", "illegal value of LDA");
        return -3;
    }
    /* Quick return */
    if (max(N, 0) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGESV, N, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetri", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB   = PLASMA_NB;

    plasma_sequence_create(plasma, &sequence);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, N, N, sequence, &request,
                             plasma_desc_mat_free(&(descA)) );
    } else {
        plasma_ziplap2tile(  descA, A, NB, NB, LDA, N, 0, 0, N, N,
                            sequence, &request);
    }

    /* Allocate workspace */
    PLASMA_Alloc_Workspace_zgetri_Tile_Async(&descA, &descW);

    /* Call the tile interface */
    PLASMA_zgetri_Tile_Async(&descA, IPIV, &descW, sequence, &request);

    if ( PLASMA_TRANSLATION == PLASMA_OUTOFPLACE ) {
        plasma_zooptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
        plasma_desc_mat_free(&descA);
    } else {
        plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
        plasma_dynamic_sync();
    }
    plasma_desc_mat_free(&(descW));

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zgetri_Tile - Computes the inverse of a matrix using the LU factorization
 *  computed by PLASMA_zgetrf.
 *  This method inverts U and then computes inv(A) by solving the system
 *  inv(A)*L = inv(U) for inv(A).
 *  Tile equivalent of PLASMA_zgetri().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the triangular factor L or U from the
 *          factorization A = P*L*U as computed by PLASMA_zgetrf.
 *          On exit, if return value = 0, the inverse of the original
 *          matrix A.
 *
 * @param[in] IPIV
 *          The pivot indices that define the permutations
 *          as returned by PLASMA_zgetrf.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, the (i,i) element of the factor U is
 *                exactly zero; The matrix is singular 
 *                and its inverse could not be computed.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetri
 * @sa PLASMA_zgetri_Tile_Async
 * @sa PLASMA_cgetri_Tile
 * @sa PLASMA_dgetri_Tile
 * @sa PLASMA_sgetri_Tile
 * @sa PLASMA_zgetrf_Tile
 *
 ******************************************************************************/
int PLASMA_zgetri_Tile(PLASMA_desc *A, int *IPIV)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descW;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetri_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);

    /* Allocate workspace */
    PLASMA_Alloc_Workspace_zgetri_Tile_Async(A, &descW);

    PLASMA_zgetri_Tile_Async(A, IPIV, &descW, sequence, &request);
    plasma_dynamic_sync();
    plasma_desc_mat_free(&(descW));

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgetri_Tile_Async - Computes the inverse of a matrix using the LU 
 *  factorization computed by PLASMA_zgetrf.
 *  This method inverts U and then computes inv(A) by solving the system
 *  inv(A)*L = inv(U) for inv(A).
 *  Non-blocking equivalent of PLASMA_zgetri_Tile().
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
 * @sa PLASMA_zgetri
 * @sa PLASMA_zgetri_Tile
 * @sa PLASMA_cgetri_Tile_Async
 * @sa PLASMA_dgetri_Tile_Async
 * @sa PLASMA_sgetri_Tile_Async
 * @sa PLASMA_zgetrf_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgetri_Tile_Async(PLASMA_desc *A, int *IPIV, PLASMA_desc *W, 
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    PLASMA_desc descW;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetri_Tile_Async", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgetri_Tile_Async", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgetri_Tile_Async", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetri_Tile_Async", "invalid A descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check descriptors for correctness */
    if (plasma_desc_check(W) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetri_Tile_Async", "invalid W descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descW = *W;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgetri_Tile_Async", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
    if (max(descA.m, 0) == 0)
        return PLASMA_SUCCESS;

    plasma_dynamic_call_5(plasma_pztrtri,
        PLASMA_enum, PlasmaUpper,
        PLASMA_enum, PlasmaNonUnit,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_9(plasma_pztrsmrv,
        PLASMA_enum, PlasmaRight,
        PLASMA_enum, PlasmaLower,
        PLASMA_enum, PlasmaNoTrans,
        PLASMA_enum, PlasmaUnit,
        PLASMA_Complex64_t, (PLASMA_Complex64_t) 1.0,
        PLASMA_desc, descA,
        PLASMA_desc, descW,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    /* No need for barrier tile2row because of previous dependencies */

    /* swap */
    plasma_dynamic_call_5(
        plasma_pzlaswpc,
        PLASMA_desc, descA,
        int *,       IPIV,
        int,         -1,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    plasma_dynamic_call_3(
        plasma_pzbarrier_row2tl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    return PLASMA_SUCCESS;
}

int PLASMA_Alloc_Workspace_zgetri_Tile_Async(PLASMA_desc *A, PLASMA_desc *W)
{
  plasma_zdesc_alloc( *W, A->mb, A->nb, A->m, A->nb, 0, 0, A->m, A->nb, 
                      plasma_desc_mat_free( W ));
  return PLASMA_SUCCESS;
}
