/**
 *
 * @file sopalin_data.h
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2024-07-05
 *
 **/
#ifndef _sopalin_data_h_
#define _sopalin_data_h_

#include "models.h"

struct sopalin_data_s {
    SolverMatrix   *solvmtx;
    double        (*cpu_coefs)[PastixKernelLvl1Nbr][8];
    double        (*gpu_coefs)[PastixKernelLvl1Nbr][8];
    pastix_model_t *cpu_models;
    pastix_model_t *gpu_models;
};
typedef struct sopalin_data_s sopalin_data_t;

void sopalin_ztrsm( pastix_data_t *pastix_data, pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_ctrsm( pastix_data_t *pastix_data, pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_dtrsm( pastix_data_t *pastix_data, pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_strsm( pastix_data_t *pastix_data, pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );

void sopalin_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_cdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_ddiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );
void sopalin_sdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data, pastix_rhs_t rhsb );

void sopalin_zgetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_cgetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_dgetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_sgetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );

void sopalin_zhetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_chetrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );

void sopalin_zpotrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_cpotrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_dpotrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_spotrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );

void sopalin_zpxtrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_cpxtrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );

void sopalin_zsytrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_csytrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_dsytrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );
void sopalin_ssytrf( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data );

#endif /* _sopalin_data_h_ */


