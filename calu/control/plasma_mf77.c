/**
 *
 * @file plasma_mf77.c
 *
 *  PLASMA mixed-precision computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Bilel Hadri
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

#define PLASMA_ZCGESV           PLASMA_FNAME(zcgesv,   ZCGESV)
#define PLASMA_DSGESV           PLASMA_FNAME(dsgesv,   DSGESV)
#define PLASMA_ZCPOSV           PLASMA_FNAME(zcposv,   ZCPOSV)
#define PLASMA_DSPOSV           PLASMA_FNAME(dsposv,   DSPOSV)
#define PLASMA_ZCGELS           PLASMA_FNAME(zcgels,   ZCGELS)
#define PLASMA_DSGELS           PLASMA_FNAME(dsgels,   DSGELS)
#define PLASMA_ZCUNGESV         PLASMA_FNAME(zcungesv, ZCUNGESV)
#define PLASMA_DSUNGESV         PLASMA_FNAME(dsungesv, DSUNGESV)

#define PLASMA_ZCGESV_TILE       PLASMA_TILE_FNAME(zcgesv,   ZCGESV)
#define PLASMA_DSGESV_TILE       PLASMA_TILE_FNAME(dsgesv,   DSGESV)
#define PLASMA_ZCPOSV_TILE       PLASMA_TILE_FNAME(zcposv,   ZCPOSV)
#define PLASMA_DSPOSV_TILE       PLASMA_TILE_FNAME(dsposv,   DSPOSV)
#define PLASMA_ZCGELS_TILE       PLASMA_TILE_FNAME(zcgels,   ZCGELS)
#define PLASMA_DSGELS_TILE       PLASMA_TILE_FNAME(dsgels,   DSGELS)
#define PLASMA_ZCUNGESV_TILE     PLASMA_TILE_FNAME(zcungesv, ZCUNGESV)
#define PLASMA_DSUNGESV_TILE     PLASMA_TILE_FNAME(dsungesv, DSUNGESV)

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (simple interface)
 **/
void PLASMA_ZCGESV(int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, int *IPIV, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_zcgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }

void PLASMA_DSGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_dsgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }

void PLASMA_ZCPOSV(PLASMA_enum *uplo, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_zcposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

void PLASMA_DSPOSV(PLASMA_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_dsposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

/* void PLASMA_ZCGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t **T, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *X, int *LDX, int *ITER, int *INFO) */
/* {   *INFO = PLASMA_zcgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); } */

/* void PLASMA_DSGELS(PLASMA_enum *trans, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO) */
/* {   *INFO = PLASMA_dsgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); } */

void PLASMA_ZCUNGESV(PLASMA_enum *trans, int *N, int *NRHS, PLASMA_Complex64_t *A, int *LDA, PLASMA_Complex64_t *B, int *LDB, PLASMA_Complex64_t *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_zcungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

void PLASMA_DSUNGESV(PLASMA_enum *trans, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
{   *INFO = PLASMA_dsungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

#ifdef __cplusplus
}
#endif

/***************************************************************************//**
 *  FORTRAN API - math functions (native interface)
 **/
void PLASMA_ZCGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_zcgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }

void PLASMA_DSGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_zcgesv_Tile((PLASMA_desc *)(*A), IPIV, (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }

void PLASMA_ZCPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_zcposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }

void PLASMA_DSPOSV_TILE(PLASMA_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_dsposv_Tile(*uplo, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }

/* void PLASMA_ZCGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO) */
/* {   *INFO = PLASMA_zcgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T), (PLASMA_desc *)(*X), ITER); } */

/* void PLASMA_DSGELS_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO) */
/* {   *INFO = PLASMA_dsgels_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*B), (PLASMA_desc *)(*T), (PLASMA_desc *)(*X), ITER); } */

void PLASMA_ZCUNGESV_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_zcungesv_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }

void PLASMA_DSUNGESV_TILE(PLASMA_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
{   *INFO = PLASMA_dsungesv_Tile(*trans, (PLASMA_desc *)(*A), (PLASMA_desc *)(*T), (PLASMA_desc *)(*B), (PLASMA_desc *)(*X), ITER); }
