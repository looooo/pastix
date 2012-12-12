/**
 *
 * @file zplrnt.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
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
 *  PLASMA_zplrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_zplrnt_Tile
 * @sa PLASMA_zplrnt_Tile_Async
 * @sa PLASMA_cplrnt
 * @sa PLASMA_dplrnt
 * @sa PLASMA_splrnt
 * @sa PLASMA_zplghe
 * @sa PLASMA_zplgsy
 *
 ******************************************************************************/
int PLASMA_zplrnt( int M, int N,
                   PLASMA_Complex64_t *A, int LDA,
                   unsigned long long int seed )
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zplrnt", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_zplrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_zplrnt", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_zplrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_ZGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zplrnt", "plasma_tune() failed");
        return status;
    }
    
    /* Set NT */
    NB = PLASMA_NB;
    plasma_sequence_create(plasma, &sequence);
    descA = plasma_desc_init(
        PlasmaComplexDouble, NB, NB, NB*NB,
        LDA, N, 0, 0, M, N);
    descA.mat = A;

    /* Call the tile interface */
    PLASMA_zplrnt_Tile_Async( &descA, seed, sequence, &request );

    plasma_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
    plasma_dynamic_sync();

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile
 *
 *  PLASMA_zplrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of PLASMA_zplrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_zplrnt
 * @sa PLASMA_zplrnt_Tile_Async
 * @sa PLASMA_cplrnt_Tile
 * @sa PLASMA_dplrnt_Tile
 * @sa PLASMA_splrnt_Tile
 * @sa PLASMA_zplghe_Tile
 * @sa PLASMA_zplgsy_Tile
 *
 ******************************************************************************/
int PLASMA_zplrnt_Tile( PLASMA_desc *A,
                        unsigned long long int seed )
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zplrnt_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_zplrnt_Tile_Async( A, seed, sequence, &request );
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex64_t_Tile_Async
 *
 *  PLASMA_zplrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of PLASMA_zplrnt_Tile().
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
 * @sa PLASMA_zplrnt
 * @sa PLASMA_zplrnt_Tile
 * @sa PLASMA_cplrnt_Tile_Async
 * @sa PLASMA_dplrnt_Tile_Async
 * @sa PLASMA_splrnt_Tile_Async
 * @sa PLASMA_zplghe_Tile_Async
 * @sa PLASMA_zplgsy_Tile_Async
 *
 ******************************************************************************/
int PLASMA_zplrnt_Tile_Async( PLASMA_desc     *A,
                              unsigned long long int seed,
                              PLASMA_sequence *sequence, 
                              PLASMA_request  *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_zplrnt_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_zplrnt_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_zplrnt_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_zplrnt_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_zplrnt_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (min( descA.m, descA.n ) == 0)
        return PLASMA_SUCCESS;

    plasma_parallel_call_4(plasma_pzplrnt,
        PLASMA_desc,      descA,
        unsigned long long int, seed,
        PLASMA_sequence*, sequence,
        PLASMA_request*,  request);

    return PLASMA_SUCCESS;
}
