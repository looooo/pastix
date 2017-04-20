/**
 * @file pastix_zcores.h
 *
 * Copyright (c) 2016      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define pastix_cblk_lock( cblk_ )    pastix_atomic_lock( &((cblk_)->lock) )
#define pastix_cblk_unlock( cblk_ )  pastix_atomic_unlock( &((cblk_)->lock) )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @addtogroup kernel_blas_lapack
 * @{
 *    This module contains all the BLAS and LAPACK-like kernels that are working
 *    on lapack layout matrices.
 */
void core_zplrnt( int m, int n, pastix_complex64_t *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );
void core_zgetro( int m, int n,
                  const pastix_complex64_t *A, int lda,
                  pastix_complex64_t *B, int ldb );
int  core_zgeadd( pastix_trans_t trans, pastix_int_t M, pastix_int_t N,
                  pastix_complex64_t alpha, const pastix_complex64_t *A, pastix_int_t LDA,
                  pastix_complex64_t beta,        pastix_complex64_t *B, pastix_int_t LDB );
int  core_zgemdm( pastix_trans_t transA, pastix_trans_t transB, int M, int N, int K,
                  pastix_complex64_t  alpha, const pastix_complex64_t *A, int LDA,
                  const pastix_complex64_t *B, int LDB,
                  pastix_complex64_t  beta, pastix_complex64_t *C, int LDC,
                  const pastix_complex64_t *D, int incD,
                  pastix_complex64_t *WORK, int LWORK );
int  core_zrrqr ( pastix_int_t m, pastix_int_t n,
                  pastix_complex64_t *A, pastix_int_t lda,
                  pastix_int_t *jpvt, pastix_complex64_t *tau,
                  pastix_complex64_t *work, pastix_int_t ldwork, double *rwork,
                  double tol, pastix_int_t nb, pastix_int_t maxrank);

void core_zpotrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zgetrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria );
void core_zhetrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria, pastix_complex64_t *work );
void core_zsytrfsp( pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t *nbpivot, double criteria, pastix_complex64_t *work );

/**
 * @}
 *
 * @addtogroup kernel_lr
 * @{
 *    This module contains all the low-rank kernels working on pastix_lr_t
 *    matrix representations.
 */
void core_zlralloc( pastix_int_t M, pastix_int_t N, pastix_int_t rkmax, pastix_lrblock_t *A );
void core_zlrfree ( pastix_lrblock_t *A );
int  core_zlrsze  ( int copy, pastix_int_t M, pastix_int_t N, pastix_lrblock_t *A, int newrk, int newrkmax );
int  core_zlr2ge  ( pastix_int_t M, pastix_int_t N, const pastix_lrblock_t *Alr, pastix_complex64_t *A, pastix_int_t lda );
int  core_zgradd  ( const pastix_lr_t *lowrank, pastix_complex64_t alpha,
                    pastix_int_t M1, pastix_int_t N1, const pastix_complex64_t *A, pastix_int_t lda,
                    pastix_int_t M2, pastix_int_t N2, pastix_lrblock_t *B,
                    pastix_int_t offx, pastix_int_t offy );
void core_zlrmm   ( const pastix_lr_t *lowrank, pastix_trans_t transA, pastix_trans_t transB,
                    pastix_int_t M, pastix_int_t N, pastix_int_t K, pastix_int_t Cm, pastix_int_t Cn,
                    pastix_int_t offx, pastix_int_t offy,
                    pastix_complex64_t alpha, const pastix_lrblock_t *A, const pastix_lrblock_t *B,
                    pastix_complex64_t beta,  pastix_lrblock_t *C,
                    pastix_complex64_t *work, pastix_int_t ldwork, SolverCblk *fcblk );
void core_zlrmge  ( const pastix_lr_t *lowrank, pastix_trans_t transA, pastix_trans_t transB,
                    pastix_int_t M, pastix_int_t N, pastix_int_t K,
                    pastix_complex64_t alpha, const pastix_lrblock_t *A, const pastix_lrblock_t *B,
                    pastix_complex64_t beta,  pastix_complex64_t *C, int ldc,
                    pastix_complex64_t *work, pastix_int_t ldwork, SolverCblk *fcblk );

/**
 *@}
 *
 * @addtogroup kernel_lr_svd
 * @{
 *    This is the SVD implementation of the low-rank kernels based on the LAPACK
 *    GESVD function.
 */
void core_zge2lr_SVD( double tol, pastix_int_t m, pastix_int_t n,
                      const pastix_complex64_t *A, pastix_int_t lda, pastix_lrblock_t *Alr );
int  core_zrradd_SVD( double tol, pastix_trans_t transA1, pastix_complex64_t alpha,
                      pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                      pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                      pastix_int_t offx, pastix_int_t offy );

void core_zge2lr_SVD_interface( pastix_fixdbl_t tol, pastix_int_t m, pastix_int_t n,
                                const void *Aptr, pastix_int_t lda, void *Alr );
int  core_zrradd_SVD_interface( pastix_fixdbl_t tol, pastix_trans_t transA1, const void *alphaptr,
                                pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                pastix_int_t offx, pastix_int_t offy );

/**
 *@}
 *
 * @addtogroup kernel_lr_rrqr
 * @{
 *    This is the rank-revealing QR implementation of the low-rank kernels based
 *    on the modified LAPACK GEQP3 function.
 */
void core_zge2lr_RRQR( double tol, pastix_int_t m, pastix_int_t n,
                       const pastix_complex64_t *A, pastix_int_t lda, pastix_lrblock_t *Alr );
int  core_zrradd_RRQR( double tol, pastix_trans_t transA1, pastix_complex64_t alpha,
                       pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                       pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                       pastix_int_t offx, pastix_int_t offy );


void core_zge2lr_RRQR_interface( pastix_fixdbl_t tol, pastix_int_t m, pastix_int_t n,
                                 const void *Aptr, pastix_int_t lda, void *Alr );

int  core_zrradd_RRQR_interface( pastix_fixdbl_t tol, pastix_trans_t transA1, const void *alphaptr,
                                 pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                 pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                 pastix_int_t offx, pastix_int_t offy );
/**
 * @}
 *
 * @addtogroup kernel_fact
 * @{
 *    This module contains all the kernel working at the solver matrix structure
 *    level for the numerical factorization step.
 */

int  core_zgeaddsp1d( const SolverCblk *cblk1, SolverCblk *cblk2,
                      const pastix_complex64_t *L1, pastix_complex64_t *L2,
                      const pastix_complex64_t *U1, pastix_complex64_t *U2 );

void cpucblk_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                      const pastix_complex64_t *A, const pastix_complex64_t *B, pastix_complex64_t *C,
                      pastix_complex64_t *work, const pastix_lr_t *lowrank );
void cpucblk_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                      pastix_trans_t trans, pastix_diag_t diag, SolverCblk *cblk,
                      const pastix_complex64_t *A, pastix_complex64_t *C, const pastix_lr_t *lowrank );

void cpublok_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, SolverCblk *fcblk,
                      pastix_int_t blok_mk, pastix_int_t blok_nk, pastix_int_t blok_mn,
                      const pastix_complex64_t *A, const pastix_complex64_t *B, pastix_complex64_t *C,
                      const pastix_lr_t *lowrank );
void cpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                      pastix_trans_t trans, pastix_diag_t diag,
                      SolverCblk *cblk, pastix_int_t blok_m,
                      const pastix_complex64_t *A, pastix_complex64_t *C,
                      const pastix_lr_t *lowrank );

#if defined(PASTIX_WITH_CUDA)
#include <cuda.h>
#include <cuComplex.h>

void gpucblk_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, const SolverBlok *blok, SolverCblk *fcblk,
                      const cuDoubleComplex *A, const cuDoubleComplex *B, cuDoubleComplex *C,
                      const pastix_lr_t *lowrank, cudaStream_t stream );

void gpublok_zgemmsp( pastix_coefside_t sideA, pastix_coefside_t sideB, pastix_trans_t trans,
                      const SolverCblk *cblk, SolverCblk *fcblk,
                      pastix_int_t blok_mk, pastix_int_t blok_nk, pastix_int_t blok_mn,
                      const cuDoubleComplex *A, const cuDoubleComplex *B, cuDoubleComplex *C,
                      const pastix_lr_t *lowrank, cudaStream_t stream );

void gpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                      pastix_trans_t trans, pastix_diag_t diag,
                      SolverCblk *cblk, pastix_int_t blok_m,
                      const cuDoubleComplex *A, cuDoubleComplex *C,
                      const pastix_lr_t *lowrank, cudaStream_t stream );

void gpu_zgemmsp_fermi( const SolverMatrix *solvmatr,
                        pastix_uplo_t uplo, pastix_trans_t trans,
                        int *blocktab,
                        const SolverCblk      *cblk,
                        const SolverBlok      *blok,
                        SolverCblk      *fcblk,
                        const cuDoubleComplex *A,
                        const cuDoubleComplex *B,
                        cuDoubleComplex *C,
                        cudaStream_t stream );
#endif

int core_zgetrfsp1d_getrf( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *U, double criteria );
int core_zgetrfsp1d_panel( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *U, double criteria,
                           const pastix_lr_t *lowrank );
int core_zgetrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk, double criteria, pastix_complex64_t *work );

int core_zpotrfsp1d_potrf( SolverCblk *cblk, pastix_complex64_t *L, double criteria );
int core_zpotrfsp1d_panel( SolverCblk *cblk, pastix_complex64_t *L, double criteria,
                           const pastix_lr_t *lowrank );
int core_zpotrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk, double criteria, pastix_complex64_t *work );

#if defined(PRECISION_z) || defined(PRECISION_c)
int core_zhetrfsp1d_hetrf( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *DLh, double criteria );
int core_zhetrfsp1d_panel( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *DLh, double criteria,
                           const pastix_lr_t *lowrank );
int core_zhetrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk, double criteria,
                           pastix_complex64_t *work1, pastix_complex64_t *work2 );

int  core_zhetrfsp1d_trsm( SolverCblk *cblk, pastix_complex64_t *L );
void core_zhetrfsp1d_gemm( SolverCblk *cblk, SolverBlok *blok, SolverCblk *fcblk,
                           pastix_complex64_t *L, pastix_complex64_t *C,
                           pastix_complex64_t *work1, pastix_complex64_t *work2 );
#endif /*defined(PRECISION_z) || defined(PRECISION_c) */

int core_zsytrfsp1d_sytrf( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *DLt, double criteria );
int core_zsytrfsp1d_panel( SolverCblk *cblk, pastix_complex64_t *L, pastix_complex64_t *DLt, double criteria,
                           const pastix_lr_t *lowrank );
int core_zsytrfsp1d      ( SolverMatrix *solvmtx, SolverCblk *cblk, double criteria,
                           pastix_complex64_t *work1, pastix_complex64_t *work2 );

int  core_zsytrfsp1d_trsm( SolverCblk *cblk, pastix_complex64_t *L );
void core_zsytrfsp1d_gemm( SolverCblk *cblk, SolverBlok *blok, SolverCblk *fcblk,
                           pastix_complex64_t *L, pastix_complex64_t *C,
                           pastix_complex64_t *work1, pastix_complex64_t *work2 );

/**
 * @}
 *
 * @addtogroup kernel_solve
 * @{
 *    This module contains all the kernel working on the solver matrix structure
 *    for the solve step.
 */

void solve_ztrsmsp( pastix_side_t side, pastix_uplo_t uplo, pastix_trans_t trans, pastix_diag_t diag,
                    SolverMatrix *datacode, SolverCblk *cblk,
                    int nrhs, pastix_complex64_t *b, int ldb );

/**
 * @}
 */
#endif /* _CORE_Z_H_ */
