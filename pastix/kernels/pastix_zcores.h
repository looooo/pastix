/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#include "z_solver.h"

#define pastix_cblk_lock( cblk_ )
#define pastix_cblk_unlock( cblk_ )
/* REQUIRED ? */
#if (!defined PRECISION_z) && (!defined PRECISION_c) && (!defined PRECISION_d) && (!defined PRECISION_s)
#define PRECISION_z
#endif
void core_zgetro(int m, int n,
                 pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB);
int core_zgeaddsp1d( z_SolverCblk * cblk1,
                     z_SolverCblk * cblk2,
                     pastix_complex64_t * L1,
                     pastix_complex64_t * L2,
                     pastix_complex64_t * U1,
                     pastix_complex64_t * U2 );

int core_zgemdm(int transA, int transB,
                int M, int N, int K,
                pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                const pastix_complex64_t *B, int LDB,
                pastix_complex64_t beta,
                pastix_complex64_t *C, int LDC,
                const pastix_complex64_t *D, int incD,
                pastix_complex64_t *WORK, int LWORK);


int core_zgetrfsp1d_getrf( z_SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_trsm( z_SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U);

int core_zgetrfsp1d( z_SolverCblk         *cblk,
                     pastix_complex64_t *L,
                     pastix_complex64_t *U,
                     double              criteria);

void core_zgetrfsp1d_gemm( z_SolverCblk         *cblk,
                           z_SolverBlok         *blok,
                           z_SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           pastix_complex64_t *Cl,
                           pastix_complex64_t *Cu,
                           pastix_complex64_t *work );

#ifdef PRECISION_z
int core_zhetrfsp1d_hetrf( z_SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work );

int core_zhetrfsp1d_trsm( z_SolverCblk         *cblk,
                          pastix_complex64_t *L);

int core_zhetrfsp1d( z_SolverCblk         *cblk,
                     pastix_complex64_t *L,
                     double              criteria,
                     pastix_complex64_t *work );

void core_zhetrfsp1d_gemm( z_SolverCblk         *cblk,
                           z_SolverBlok         *blok,
                           z_SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2 );
#endif /* PRECISION_z */
int core_zpotrfsp1d_potrf( z_SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria );

int core_zpotrfsp1d_trsm( z_SolverCblk         *cblk,
                          pastix_complex64_t *L );

int core_zpotrfsp1d( z_SolverCblk         *cblk,
                     pastix_complex64_t *L,
                     double              criteria );

void core_zpotrfsp1d_gemm(z_SolverCblk         *cblk,
                          z_SolverBlok         *blok,
                          z_SolverCblk         *fcblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work);

int core_zsytrfsp1d_sytrf( z_SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work );

int core_zsytrfsp1d_trsm( z_SolverCblk         *cblk,
                          pastix_complex64_t *L );

int core_zsytrfsp1d( z_SolverCblk         *cblk,
                     pastix_complex64_t *L,
                     double              criteria,
                     pastix_complex64_t *work );

void core_zsytrfsp1d_gemm( z_SolverCblk         *cblk,
                           z_SolverBlok         *blok,
                           z_SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2);

#endif /* _CORE_Z_H_ */
