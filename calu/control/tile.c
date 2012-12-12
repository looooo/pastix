/**
 *
 * @file tile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include "common.h"
#include "auxiliary.h"
#include "tile.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Lapack_to_Tile(void *Af77, int LDA, PLASMA_desc *A)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Lapack_to_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_Lapack_to_Tile", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    switch( A->dtyp ) {
    case PlasmaRealFloat:
      plasma_parallel_call_5(
          plasma_pslapack_to_tile,
          float*, Af77,
          int, LDA,
          PLASMA_desc, descA,
          PLASMA_sequence*, sequence,
          PLASMA_request*, &request);
      break;
    case PlasmaRealDouble:
        plasma_parallel_call_5(
            plasma_pdlapack_to_tile,
            double*, Af77,
            int, LDA,
            PLASMA_desc, descA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
      break;
    case PlasmaComplexFloat:
        plasma_parallel_call_5(
            plasma_pclapack_to_tile,
            PLASMA_Complex32_t*, Af77,
            int, LDA,
            PLASMA_desc, descA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
        break;
    case PlasmaComplexDouble:
        plasma_parallel_call_5(
            plasma_pzlapack_to_tile,
            PLASMA_Complex64_t*, Af77,
            int, LDA,
            PLASMA_desc, descA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
        break;
    default:
        plasma_error("PLASMA_Lapack_to_Tile", "Type unknown");
    }
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the PLASMA matrix in tile layout.
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Tile_to_Lapack(PLASMA_desc *A, void *Af77, int LDA)
{
    PLASMA_desc descA = *A;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Tile_to_Lapack", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (plasma_desc_check(&descA) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_Tile_to_Lapack", "invalid descriptor");
        return PLASMA_ERR_ILLEGAL_VALUE;
    }
    plasma_sequence_create(plasma, &sequence);
    switch( A->dtyp ) {
    case PlasmaRealFloat:
      plasma_parallel_call_5(
          plasma_pstile_to_lapack,
          PLASMA_desc, descA,
          float*, Af77,
          int, LDA,
          PLASMA_sequence*, sequence,
          PLASMA_request*, &request);
      break;
    case PlasmaRealDouble:
        plasma_parallel_call_5(
            plasma_pdtile_to_lapack,
            PLASMA_desc, descA,
            double*, Af77,
            int, LDA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
      break;
    case PlasmaComplexFloat:
        plasma_parallel_call_5(
            plasma_pctile_to_lapack,
            PLASMA_desc, descA,
            PLASMA_Complex32_t*, Af77,
            int, LDA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
        break;
    case PlasmaComplexDouble:
        plasma_parallel_call_5(
            plasma_pztile_to_lapack,
            PLASMA_desc, descA,
            PLASMA_Complex64_t*, Af77,
            int, LDA,
            PLASMA_sequence*, sequence,
            PLASMA_request*, &request);
        break;
    default:
        plasma_error("PLASMA_Tile_to_Lapack", "Type unknown");
    }
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}
