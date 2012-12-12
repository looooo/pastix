/**
 *
 * @file zgetrf.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2009-11-15
 *
 * @precisions normal z -> s d c
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t
 *
 *  PLASMA_zgetrf - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm with partial tile pivoting with row interchanges.
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
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrf_Tile
 * @sa PLASMA_zgetrf_Tile_Async
 * @sa PLASMA_cgetrf
 * @sa PLASMA_dgetrf
 * @sa PLASMA_sgetrf
 *
 ******************************************************************************/
int PLASMA_zgetrf(int M, int N,
                  PLASMA_Complex64_t *A, int LDA,
                  int *IPIV)
{
    int NB, NBNB, minMN;
    int status;
    PLASMA_desc descA ;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zgetrf", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zgetrf", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zgetrf", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = plasma_tune(PLASMA_FUNC_ZGESV, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrf", "plasma_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB   = PLASMA_NB;
    NBNB = NB*NB;

    plasma_sequence_create(plasma, &sequence);

    descA = plasma_desc_init(
        PlasmaComplexDouble,
        NB, NB, NBNB,
        LDA, N, 0, 0, M, N);
    descA.mat = A;

    minMN = min(M, N);
    memset(IPIV, 0, minMN*sizeof(int));

    /* Call the tile interface */
    plasma_dynamic_call_4(plasma_pzgetrf_reclap,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, &request);

    plasma_dynamic_sync();

    /*
     * Generate the correct IPIV (Has to be moved in a task)
     */
    {
        int i, inc, tmp, j;
        for(i=1; i<descA.mt; i++) {
            inc = i*descA.mb;
            tmp = min( minMN - inc, descA.mb);
            if ( tmp < 1 )
              break;

            for (j=0; j<tmp; j++)
                IPIV[inc+j] = IPIV[inc+j] + inc;
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
 *  PLASMA_zgetrf_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of PLASMA_zgetrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[out] IPIV
 *          The pivot indices that define the permutations
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, and division by zero will occur
 *               if it is used to solve a system of equations.
 *
 *******************************************************************************
 *
 * @sa PLASMA_zgetrf
 * @sa PLASMA_zgetrf_Tile_Async
 * @sa PLASMA_cgetrf_Tile
 * @sa PLASMA_dgetrf_Tile
 * @sa PLASMA_sgetrf_Tile
 * @sa PLASMA_zgetrs_Tile
 *
 ******************************************************************************/
int PLASMA_zgetrf_Tile(PLASMA_desc *A, int *IPIV)
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zgetrf_Tile_Async(A, IPIV, sequence, &request);
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zgetrf_Tile_Async - Computes the tile LU factorization of a matrix.
 *  Non-blocking equivalent of PLASMA_zgetrf_Tile().
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
 * @sa PLASMA_zgetrf
 * @sa PLASMA_zgetrf_Tile
 * @sa PLASMA_cgetrf_Tile_Async
 * @sa PLASMA_dgetrf_Tile_Async
 * @sa PLASMA_sgetrf_Tile_Async
 * @sa PLASMA_zgetrs_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zgetrf_Tile_Async(PLASMA_desc *A, int *IPIV,
                             PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zgetrf_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zgetrf_Tile", "invalid first descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zgetrf_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;
*/

#if 1
    PLASMA_desc descW;
    int Wm = ((descA.mt + 3) / 4 ) * 4 * descA.mb;
    int tutu = Wm, toto = 1;
    while ( (tutu >> 2) > 1 ){
        toto++;
        tutu = tutu >> 2;
    }
    plasma_zdesc_alloc( descW, 4*descA.mb, descA.nb, Wm, toto*descA.nb, 0, 0, Wm, toto*descA.nb, (void)descW);

    /* int *Wi = (int *)malloc(Wm * 2 * sizeof(int)); */
    /* memset(Wi, 0, Wm * 2 * sizeof(int)); */
    int *Wi = (int *)malloc(Wm * toto * sizeof(int));
    memset(Wi, 0, Wm * toto * sizeof(int));

    plasma_dynamic_call_6(
        plasma_pzgetrf_tntpiv,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_desc, descW,
        int*, Wi,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_sync();
    plasma_desc_mat_free(&descW);
    free(Wi);
#else
    plasma_dynamic_call_3(
        plasma_pzbarrier_tl2pnl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    plasma_dynamic_call_4(plasma_pzgetrf_rectil,
        PLASMA_desc, descA,
        int*, IPIV,
        PLASMA_sequence*, sequence,
        PLASMA_request*, request);

    plasma_dynamic_call_3(
        plasma_pzbarrier_pnl2tl,
        PLASMA_desc, descA,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);
#endif

    return PLASMA_SUCCESS;
}
