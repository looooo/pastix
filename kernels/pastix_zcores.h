/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#define pastix_cblk_lock( cblk_ )
#define pastix_cblk_unlock( cblk_ )

pastix_int_t core_z_compress_LR(pastix_complex64_t *fL,
                                pastix_int_t stride,
                                pastix_int_t dimb,
                                pastix_int_t dima,
                                pastix_complex64_t *u,
                                pastix_int_t ldu,
                                pastix_complex64_t *v,
                                pastix_int_t ldv);

void core_z_uncompress_LR(pastix_complex64_t *fL,
                          pastix_int_t stride,
                          pastix_int_t dimb,
                          pastix_int_t dima,
                          pastix_complex64_t *u,
                          pastix_int_t ldu,
                          pastix_complex64_t *v,
                          pastix_int_t ldv,
                          pastix_int_t rank);

pastix_int_t core_z_add_LR(pastix_complex64_t *u1,
                           pastix_complex64_t *v1,
                           pastix_int_t dim_u1,
                           pastix_int_t dim_v1,
                           pastix_int_t rank_1,
                           pastix_int_t ld_u1,
                           pastix_int_t ld_v1,
                           pastix_complex64_t *u2,
                           pastix_complex64_t *v2,
                           pastix_int_t dim_u2,
                           pastix_int_t dim_v2,
                           pastix_int_t rank_2,
                           pastix_int_t ld_u2,
                           pastix_int_t ld_v2,
                           pastix_int_t x2,
                           pastix_int_t y2);

void core_zgetro(int m, int n,
                 pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB);

int core_zgeaddsp1d( SolverCblk * cblk1,
                     SolverCblk * cblk2,
                     pastix_complex64_t * L1,
                     pastix_complex64_t * L2,
                     pastix_complex64_t * U1,
                     pastix_complex64_t * U2 );

int core_zgemdm(int transA, int transB,
                int M, int N, int K,
                      pastix_complex64_t  alpha,
                const pastix_complex64_t *A,    int LDA,
                const pastix_complex64_t *B,    int LDB,
                      pastix_complex64_t  beta,
                      pastix_complex64_t *C,    int LDC,
                const pastix_complex64_t *D,    int incD,
                      pastix_complex64_t *WORK, int LWORK);

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U);

/* TODO: add properly */
int core_zgetrfsp1d_trsm2( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U);

int core_zgetrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

void core_zgetrfsp1d_gemm( SolverCblk         *cblk,
                           SolverBlok         *blok,
                           SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           pastix_complex64_t *Cl,
                           pastix_complex64_t *Cu,
                           pastix_complex64_t *work );

int core_zgetrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work );

#if defined(PRECISION_z) || defined(PRECISION_c)
int core_zhetrfsp1d_hetrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work );

int core_zhetrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L);

int core_zhetrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work);

void core_zhetrfsp1d_gemm( SolverCblk         *cblk,
                           SolverBlok         *blok,
                           SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2 );

int core_zhetrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work1,
                     pastix_complex64_t *work2 );

#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

int core_zpotrfsp1d_potrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria );

int core_zpotrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L );

int core_zpotrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria);

void core_zpotrfsp1d_gemm(SolverCblk         *cblk,
                          SolverBlok         *blok,
                          SolverCblk         *fcblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work);

int core_zpotrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work );

int core_zsytrfsp1d_sytrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work );

int core_zsytrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L );

int core_zsytrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work);

void core_zsytrfsp1d_gemm( SolverCblk         *cblk,
                           SolverBlok         *blok,
                           SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2);

int core_zsytrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work1,
                     pastix_complex64_t *work2 );

#endif /* _CORE_Z_H_ */
