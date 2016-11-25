/**
 * @file pastix_zcores.h
 *
 * Copyright (c) 2016      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#define pastix_cblk_lock( cblk_ )    pastix_atomic_lock( &((cblk_)->lock) )
#define pastix_cblk_unlock( cblk_ )  pastix_atomic_unlock( &((cblk_)->lock) )

int
core_zlralloc( pastix_int_t M, pastix_int_t N,
               pastix_int_t rkmax, pastix_lrblock_t *A );

int
core_zlrfree( pastix_lrblock_t *A );

int
core_zlrsze( int copy, pastix_int_t M, pastix_int_t N,
             pastix_lrblock_t *A, int newrk, int newrkmax );

int
core_zge2lr_SVD( double tol, pastix_int_t M, pastix_int_t N,
                 const pastix_complex64_t *A, pastix_int_t lda,
                 pastix_lrblock_t *Alr );

int
core_zlr2ge( pastix_int_t M, pastix_int_t N,
             const pastix_lrblock_t *Alr,
             pastix_complex64_t *A, pastix_int_t lda );

int
core_zrradd_SVD( double tol, int transA1, pastix_complex64_t alpha,
                 pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                 pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                 pastix_int_t offx, pastix_int_t offy );

int
core_zrradd_RRQR( double tol, int transA1, pastix_complex64_t alpha,
                  pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                  pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                  pastix_int_t offx, pastix_int_t offy );

int
core_zrradd( double tol, int transA1, pastix_complex64_t alpha,
             pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
             pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
             pastix_int_t offx, pastix_int_t offy );

int
core_zgradd( double tol, pastix_complex64_t alpha,
             pastix_int_t M1, pastix_int_t N1, pastix_complex64_t *A, pastix_int_t lda,
             pastix_int_t M2, pastix_int_t N2, pastix_lrblock_t   *B,
             pastix_int_t offx, pastix_int_t offy );

int
core_zlrmm( double tol, int transA, int transB,
            pastix_int_t M, pastix_int_t N, pastix_int_t K,
            pastix_int_t Cm, pastix_int_t Cn,
            pastix_int_t offx, pastix_int_t offy,
            pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                      const pastix_lrblock_t *B,
            pastix_complex64_t beta,  pastix_lrblock_t *C,
            pastix_complex64_t *work, pastix_int_t ldwork,
            SolverCblk *fcblk );

int
core_zlrmge( double tol, int transA, int transB,
             pastix_int_t M, pastix_int_t N, pastix_int_t K,
             pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                       const pastix_lrblock_t *B,
             pastix_complex64_t beta,  pastix_complex64_t *C, int ldc,
             pastix_complex64_t *work, pastix_int_t ldwork,
             SolverCblk *fcblk );

void core_zplrnt( int m, int n, pastix_complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );

void core_zgetro(int m, int n,
                 const pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

void core_zgetrox(pastix_complex64_t alpha, int m, int n,
                  const pastix_complex64_t *A, int lda,
                  pastix_complex64_t *B, int ldb);

int core_zgeadd( pastix_int_t trans, pastix_int_t M, pastix_int_t N,
                       pastix_complex64_t  alpha,
                 const pastix_complex64_t *A, pastix_int_t LDA,
                       pastix_complex64_t  beta,
                       pastix_complex64_t *B, pastix_int_t LDB);

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

void core_zgemmsp( int uplo, int trans,
                   const SolverCblk         *cblk,
                   const SolverBlok         *blok,
                         SolverCblk         *fcblk,
                   const pastix_complex64_t *A,
                   const pastix_complex64_t *B,
                         pastix_complex64_t *C,
                         pastix_complex64_t *work,
                         double              tol  );

void
core_zgemmsp_2d2dsub( int uplo, int trans,
                      pastix_int_t blok_mk,
                      pastix_int_t blok_kn,
                      pastix_int_t blok_mn,
                      const SolverCblk         *cblk,
                            SolverCblk         *fcblk,
                      const pastix_complex64_t *A,
                      const pastix_complex64_t *B,
                            pastix_complex64_t *C );

void core_zgemmsp_2dlrsub( int coef,
                           int uplo, int trans,
                           pastix_int_t blok_mk,
                           pastix_int_t blok_kn,
                           pastix_int_t blok_mn,
                     const SolverCblk         *cblk,
                           SolverCblk         *fcblk,
                           double tolerance );

void core_ztrsmsp( int coef, int side, int uplo, int trans, int diag,
                   SolverCblk         *cblk,
             const pastix_complex64_t *A,
                   pastix_complex64_t *C );

int core_ztrsmsp_2dsub( int side, int uplo, int trans, int diag,
                              SolverCblk         *cblk,
                              pastix_int_t        fcblknum,
                        const pastix_complex64_t *A,
                              pastix_complex64_t *C );

int core_ztrsmsp_2dlrsub( int coef, int side, int uplo, int trans, int diag,
                          SolverCblk   *cblk,
                          pastix_int_t  blok_m );

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria);

int core_zgetrfsp1d( SolverMatrix       *solvmtx,
                     SolverCblk         *cblk,
                     double              criteria,
                     pastix_complex64_t *work,
                     double              tol );

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

void gpu_zgemmsp( int uplo, int trans,
                  const SolverCblk      *cblk,
                  const SolverBlok      *blok,
                        SolverCblk      *fcblk,
                  const cuDoubleComplex *A,
                  const cuDoubleComplex *B,
                        cuDoubleComplex *C,
                        cudaStream_t stream );
void
gpu_zgemmsp_2d2dsub( int trans,
                     pastix_int_t blok_mk,
                     pastix_int_t blok_kn,
                     pastix_int_t blok_mn,
                     const SolverCblk      *cblk,
                           SolverCblk      *fcblk,
                     const cuDoubleComplex *A,
                     const cuDoubleComplex *B,
                           cuDoubleComplex *C,
                           cudaStream_t stream );
#endif

void solve_ztrsmsp( int side, int uplo, int trans, int diag,
                    SolverMatrix *datacode, SolverCblk *cblk,
                    int nrhs, pastix_complex64_t *b, int ldb );

int
core_zrrqr( pastix_int_t m, pastix_int_t n,
            pastix_complex64_t *A, pastix_int_t lda,
            pastix_int_t *jpvt, pastix_complex64_t *tau,
            pastix_complex64_t *work, pastix_int_t ldwork,
            double *rwork,
            double tol, pastix_int_t nb, pastix_int_t maxrank );

int
core_zge2lr_RRQR( double tol, pastix_int_t m, pastix_int_t n,
                  const pastix_complex64_t *A, pastix_int_t lda,
                  pastix_lrblock_t *Alr );

int core_zge2lr( double tol, pastix_int_t m, pastix_int_t n,
                     const pastix_complex64_t *A, pastix_int_t lda,
                     void *Alr );

#endif /* _CORE_Z_H_ */
