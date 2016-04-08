/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#define pastix_cblk_lock( cblk_ )    pastix_atomic_lock( &((cblk_)->lock) )
#define pastix_cblk_unlock( cblk_ )  pastix_atomic_unlock( &((cblk_)->lock) )

#define L_side 0
#define U_side 1

void core_zlr2ge( pastix_int_t m, pastix_int_t n,
                  pastix_lrblock_t *Alr,
                  const pastix_complex64_t *A, pastix_int_t lda );
void core_zge2lr( double tol, pastix_int_t m, pastix_int_t n,
                  const pastix_complex64_t *A, pastix_int_t lda,
                  pastix_lrblock_t *Alr );

pastix_int_t
core_z_compress_LR(double tol, pastix_int_t m, pastix_int_t n,
                   const pastix_complex64_t *A, pastix_int_t lda,
                   pastix_complex64_t *u, pastix_int_t ldu,
                   pastix_complex64_t *v, pastix_int_t ldv);

void
core_z_uncompress_LR( pastix_int_t m, pastix_int_t n, pastix_int_t rank,
                      const pastix_complex64_t *u, pastix_int_t ldu,
                      const pastix_complex64_t *v, pastix_int_t ldv,
                      pastix_complex64_t *A, pastix_int_t lda );

pastix_int_t core_z_add_LR(pastix_complex64_t *u1,
                           pastix_complex64_t *v1,
                           pastix_int_t dim_u1,
                           pastix_int_t dim_v1,
                           pastix_int_t rank_1,
                           pastix_int_t ld_u1,
                           pastix_int_t ld_v1,
                           pastix_complex64_t *u2,
                           pastix_complex64_t *v2,
                           pastix_int_t trans,
                           pastix_int_t dim_u2,
                           pastix_int_t dim_v2,
                           pastix_int_t rank_2,
                           pastix_int_t ld_u2,
                           pastix_int_t ld_v2,
                           pastix_int_t x2,
                           pastix_int_t y2);

void core_z_lr2dense(SolverBlok *blok,
                     pastix_complex64_t *A,
                     pastix_int_t stride,
                     pastix_int_t width,
                     pastix_int_t side);

void core_zproduct_lr2dense(SolverBlok *blok1,
                            pastix_complex64_t *A1,
                            pastix_int_t stride1,
                            pastix_int_t width1,
                            pastix_int_t side1,
                            SolverBlok *blok2,
                            pastix_complex64_t *A2,
                            pastix_int_t stride2,
                            pastix_int_t width2,
                            pastix_int_t side2,
                            pastix_complex64_t *work,
                            pastix_int_t ldwork);

void core_zproduct_lr2lr(SolverBlok *blok1,
                         pastix_complex64_t *A1,
                         pastix_int_t stride1,
                         pastix_int_t width1,
                         pastix_int_t side1,
                         SolverBlok *blok2,
                         pastix_complex64_t *A2,
                         pastix_int_t stride2,
                         pastix_int_t width2,
                         pastix_int_t side2,
                         SolverBlok *blok3,
                         pastix_complex64_t *A3,
                         pastix_complex64_t *A3_contrib,
                         pastix_int_t stride3,
                         pastix_int_t width3,
                         pastix_int_t side3,
                         pastix_int_t offset,
                         pastix_complex64_t *work);

void core_zgetro(int m, int n,
                 pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB);

int core_zgemdm(int transA, int transB,
                int M, int N, int K,
                      pastix_complex64_t  alpha,
                const pastix_complex64_t *A,    int LDA,
                const pastix_complex64_t *B,    int LDB,
                      pastix_complex64_t  beta,
                      pastix_complex64_t *C,    int LDC,
                const pastix_complex64_t *D,    int incD,
                      pastix_complex64_t *WORK, int LWORK);

int core_zgeaddsp1d( SolverCblk * cblk1,
                     SolverCblk * cblk2,
                     pastix_complex64_t * L1,
                     pastix_complex64_t * L2,
                     pastix_complex64_t * U1,
                     pastix_complex64_t * U2 );

void core_zgemmsp( int diag, int trans,
                   SolverCblk         *cblk,
                   SolverBlok         *blok,
                   SolverCblk         *fcblk,
                   pastix_complex64_t *A,
                   pastix_complex64_t *B,
                   pastix_complex64_t *C,
                   pastix_complex64_t *work );

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U);

int core_zgetrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work );

int core_zgetrfsp1d_LR( SolverMatrix       *solvmtx,
                        SolverCblk         *cblk,
                        double              criteria,
                        pastix_complex64_t *work);

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

#if defined(PASTIX_WITH_CUDA)
#include <cuda.h>
#include <cuComplex.h>

typedef struct gemm_params_s{
    int M[32];
    int Acoefind[32];
    void *Carray[32];
} gemm_params_t;

void pastix_zgemm_vbatched_nt(
    gemm_params_t params,
    pastix_trans_t transB,
    pastix_int_t n, pastix_int_t k,
    cuDoubleComplex alpha,
    cuDoubleComplex const * dA, pastix_int_t ldda,
    cuDoubleComplex const * dB, pastix_int_t lddb,
    cuDoubleComplex beta,
    pastix_int_t lddc,
    pastix_int_t max_m, pastix_int_t batchCount, cudaStream_t stream );

void gpu_zgemmsp( int uplo, int trans,
                  SolverCblk      *cblk,
                  SolverBlok      *blok,
                  SolverCblk      *fcblk,
                  cuDoubleComplex *A,
                  cuDoubleComplex *B,
                  cuDoubleComplex *C,
                  cudaStream_t stream );
#endif

#endif /* _CORE_Z_H_ */
