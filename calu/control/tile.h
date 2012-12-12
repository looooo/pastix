/**
 *
 * @file tile.h
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
#ifndef _PLASMA_TILE_H_
#define _PLASMA_TILE_H_

#if defined( _WIN32 ) || defined( _WIN64 )
typedef __int64 int64_t;
#else
#include <inttypes.h>
#endif

#define ELTADDR(A, type, m, n)  (type *)plasma_geteltaddr(A, m, n)
#define ELTLDD(A, k) ( ( (((k)-1)/(A).mb) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define BLKADDR(A, type, m, n)  (type *)plasma_getaddr(A, m, n)
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines - static scheduling
 **/
void plasma_pztile_to_lapack(plasma_context_t *plasma);
void plasma_pctile_to_lapack(plasma_context_t *plasma);
void plasma_pdtile_to_lapack(plasma_context_t *plasma);
void plasma_pstile_to_lapack(plasma_context_t *plasma);
void plasma_pzlapack_to_tile(plasma_context_t *plasma);
void plasma_pclapack_to_tile(plasma_context_t *plasma);
void plasma_pdlapack_to_tile(plasma_context_t *plasma);
void plasma_pslapack_to_tile(plasma_context_t *plasma);
void plasma_pztile_zero(plasma_context_t *plasma);
void plasma_pctile_zero(plasma_context_t *plasma);
void plasma_pdtile_zero(plasma_context_t *plasma);
void plasma_pstile_zero(plasma_context_t *plasma);

/***************************************************************************//**
 *  Internal routines - dynamic scheduling
 **/
void plasma_pztile_to_lapack_quark(PLASMA_desc, PLASMA_Complex64_t*, int, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctile_to_lapack_quark(PLASMA_desc, PLASMA_Complex32_t*, int, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtile_to_lapack_quark(PLASMA_desc, double*, int, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstile_to_lapack_quark(PLASMA_desc, float*, int, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlapack_to_tile_quark(PLASMA_Complex64_t*, int, PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pclapack_to_tile_quark(PLASMA_Complex32_t*, int, PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdlapack_to_tile_quark(double*, int, PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pslapack_to_tile_quark(float*, int, PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztile_zero_quark(PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pctile_zero_quark(PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pdtile_zero_quark(PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pstile_zero_quark(PLASMA_desc, PLASMA_sequence *sequence, PLASMA_request *request);

#ifdef __cplusplus
}
#endif

#endif
