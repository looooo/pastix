/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef Z_COMPUTE_DIAG_H
#define Z_COMPUTE_DIAG_H

#include "sopalin_define.h"

#define z_PASTIX_getrf       API_CALL(z_PASTIX_getrf)
#define z_PASTIX_potrf       API_CALL(z_PASTIX_potrf)
#define z_PASTIX_sytrf       API_CALL(z_PASTIX_sytrf)
#define z_PASTIX_hetrf       API_CALL(z_PASTIX_hetrf)
#define z_PASTIX_getrf_block API_CALL(z_PASTIX_getrf_block)
#define z_PASTIX_potrf_block API_CALL(z_PASTIX_potrf_block)
#define z_PASTIX_sytrf_block API_CALL(z_PASTIX_sytrf_block)
#define z_PASTIX_hetrf_block API_CALL(z_PASTIX_hetrf_block)
#define z_DimTrans           API_CALL(z_DimTrans)
#define z_factor_diag        API_CALL(z_factor_diag)

void z_PASTIX_getrf       ( pastix_complex64_t *A, pastix_int_t m, pastix_int_t n, pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_potrf       ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_sytrf       ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_hetrf       ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_getrf_block ( pastix_complex64_t *A, pastix_int_t m, pastix_int_t n, pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_potrf_block ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit );
void z_PASTIX_sytrf_block ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit, pastix_complex64_t * tmp4   );
void z_PASTIX_hetrf_block ( pastix_complex64_t *A, pastix_int_t n,        pastix_int_t lda, pastix_int_t *npvt, double crit, pastix_complex64_t * tmp4   );
void z_DimTrans           ( pastix_complex64_t *A, pastix_int_t lda, pastix_int_t size, pastix_complex64_t *B );

#endif /* Z_COMPUTE_DIAG_H */
