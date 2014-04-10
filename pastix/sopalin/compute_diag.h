#ifndef _COMPUTE_DIAG_H
#define _COMPUTE_DIAG_H

#include "redefine_functions.h"
#include "sopalin_define.h"

#define PASTIX_getrf       API_CALL(PASTIX_getrf)
#define PASTIX_potrf       API_CALL(PASTIX_potrf)
#define PASTIX_sytrf       API_CALL(PASTIX_sytrf)
#define PASTIX_hetrf       API_CALL(PASTIX_hetrf)
#define PASTIX_getrf_block API_CALL(PASTIX_getrf_block)
#define PASTIX_potrf_block API_CALL(PASTIX_potrf_block)
#define PASTIX_sytrf_block API_CALL(PASTIX_sytrf_block)
#define PASTIX_hetrf_block API_CALL(PASTIX_hetrf_block)
#define DimTrans           API_CALL(DimTrans)
#define factor_diag        API_CALL(factor_diag)

void PASTIX_getrf       ( pastix_float_t *A, pastix_int_t m, pastix_int_t n, pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_potrf       ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_sytrf       ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_hetrf       ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_getrf_block ( pastix_float_t *A, pastix_int_t m, pastix_int_t n, pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_potrf_block ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void PASTIX_sytrf_block ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit, pastix_float_t * tmp4   );
void PASTIX_hetrf_block ( pastix_float_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit, pastix_float_t * tmp4   );
void DimTrans           ( pastix_float_t *A, pastix_int_t lda, pastix_int_t size, pastix_float_t *B );

#endif _COMPUTE_DIAG_H
