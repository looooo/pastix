/**
 *
 * @file z_tests.h
 *
 * Tests functions header.
 *
 * @copyright 2018-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2024-07-05
 *
 * @precisions normal z -> z c d s
 *
 **/
#ifndef _z_tests_h_
#define _z_tests_h_

#include "pastix_lowrank.h"
#include "tests.h"

extern pastix_lr_t z_lowrank;

void z_init( pastix_data_t *pastix_data, unsigned long long seed, pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda );
int z_bcsc_spmv_time( pastix_data_t *pastix_data, const spmatrix_t *spm, pastix_int_t nrhs );

int z_bcsc_spmv_check( spm_trans_t trans, const spmatrix_t *spm, pastix_data_t *pastix_data );
int z_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc );
int z_bvec_gemv_check( pastix_data_t *pastix_data, int check, int m, int n );
int z_bvec_check( pastix_data_t *pastix_data );
int z_bvec_time( pastix_data_t *pastix_data );
int z_bvec_compare( pastix_data_t            *pastix_data,
                    pastix_int_t              m,
                    pastix_int_t              n,
                    const pastix_complex64_t *A,
                    pastix_int_t              lda,
                    const pastix_complex64_t *B,
                    pastix_int_t              ldb );
int z_bvec_applyorder_check ( pastix_data_t *pastix_data,
                              spmatrix_t    *spm,
                              pastix_int_t   nrhs );

int  z_lowrank_genmat( int mode, double tolerance, double threshold, test_matrix_t *A );
void z_lowrank_genmat_comp( const pastix_lr_t *lowrank, int mode, double threshold, test_matrix_t *A );

int z_lowrank_check_ge2lr( pastix_lr_t *lowrank, test_matrix_t *C );

int z_lowrank_check_rradd( pastix_lr_t *lowrank,
                           pastix_int_t offx, pastix_int_t offy,
                           pastix_complex64_t zalpha,
                           const test_matrix_t *A,
                           const test_matrix_t *B,
                           test_matrix_t *C );

int z_lowrank_check_lrmm( pastix_lr_t *lowrank,
                          pastix_int_t offx, pastix_int_t offy,
                          pastix_complex64_t   alpha,
                          const test_matrix_t *A,
                          const test_matrix_t *B,
                          pastix_complex64_t   beta,
                          const test_matrix_t *C );

#endif /* _z_tests_h_ */
