/**
 *
 * @file core_zgelrops_update.c
 *
 * PaStiX low-rank kernel routines
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @author Pierre Ramet
 * @date 2017-11-06
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "z_nan_check.h"
#include "kernels_trace.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t mzone = -1.0;
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

void
core_zfrfr2fr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               pastix_int_t Cm, pastix_int_t Cn,
               pastix_int_t offx, pastix_int_t offy,
               pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                         const pastix_lrblock_t *B,
               pastix_complex64_t beta,        pastix_lrblock_t *C,
               pastix_complex64_t *work, pastix_int_t lwork,
               SolverCblk *fcblk )
{
    pastix_int_t ldau, ldbu, ldcu;
    pastix_complex64_t *Cptr;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldcu = Cm;

    Cptr  = C->u;
    Cptr += ldcu * offy + offx;

    pastix_cblk_lock( fcblk );
    assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */

    /*
     * Everything is full rank we apply directly a GEMM
     */
    kernel_trace_start_lvl2( PastixKernelLvl2_FRFR2FR );
    cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                 M, N, K,
                 CBLAS_SADDR(alpha), A->u, ldau,
                 B->u, ldbu,
                 CBLAS_SADDR(beta),  Cptr, ldcu );
    kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, K ) );

    pastix_cblk_unlock( fcblk );

    (void)lowrank;
    (void)Cn;
    (void)work;
    (void)lwork;
}

void
core_zfrlr2fr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               pastix_int_t Cm, pastix_int_t Cn,
               pastix_int_t offx, pastix_int_t offy,
               pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                         const pastix_lrblock_t *B,
               pastix_complex64_t beta,        pastix_lrblock_t *C,
               pastix_complex64_t *work, pastix_int_t lwork,
               SolverCblk *fcblk )
{
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldbu, ldbv, ldcu;
    int allocated = 0;

    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( M, B->rk, K ) + FLOPS_ZGEMM( M, N, B->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( K, N, B->rk ) + FLOPS_ZGEMM( M, N, K     );

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    kernel_trace_start_lvl2( PastixKernelLvl2_FRLR2FR );
    if ( flops1 <= flops2 ) {
        if ( lwork < M * B->rk ) {
            work = malloc( M * B->rk * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  (A * Bv) * Bu^t
         */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, M );

        pastix_cblk_lock( fcblk );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, B->rk,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        kernel_trace_stop_lvl2( flops1 );
        pastix_cblk_unlock( fcblk );
    }
    else {
        if ( lwork < K * N ) {
            work = malloc( K * N * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  A * (Bu * Bv^t)^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     K, N, B->rk,
                     CBLAS_SADDR(zone),  B->u, ldbu,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, K );

        pastix_cblk_lock( fcblk );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, K,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        kernel_trace_stop_lvl2( flops2 );
        pastix_cblk_unlock( fcblk );
    }

    if ( allocated ) {
        free( work );
    }

    (void)lowrank;
    (void)Cn;
}

void
core_zlrfr2fr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               pastix_int_t Cm, pastix_int_t Cn,
               pastix_int_t offx, pastix_int_t offy,
               pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                         const pastix_lrblock_t *B,
               pastix_complex64_t beta,        pastix_lrblock_t *C,
               pastix_complex64_t *work, pastix_int_t lwork,
               SolverCblk *fcblk )
{
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldav, ldbu, ldcu;
    int allocated = 0;

    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( A->rk, N, K ) + FLOPS_ZGEMM( M, N, A->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( M, K, A->rk ) + FLOPS_ZGEMM( M, N, K     );

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    kernel_trace_start_lvl2( PastixKernelLvl2_LRFR2FR );
    if ( flops1 <= flops2 ) {
        if ( lwork < A->rk * N ) {
            work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  Au * (Av^t * B^t)
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v, ldav,
                                         B->u, ldbu,
                     CBLAS_SADDR(zzero), work, A->rk );

        pastix_cblk_lock( fcblk );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, N, A->rk,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, A->rk,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        kernel_trace_stop_lvl2( flops1 );
        pastix_cblk_unlock( fcblk );
    }
    else {
        if ( lwork < M * K ) {
            work = malloc( M * K * sizeof(pastix_complex64_t) );
            allocated = 1;
        }

        /*
         *  (Au * Av^t) * B^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, K, A->rk,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         A->v, ldav,
                     CBLAS_SADDR(zzero), work, M );

        pastix_cblk_lock( fcblk );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        kernel_trace_stop_lvl2( flops2 );
        pastix_cblk_unlock( fcblk );
    }

    if ( allocated ) {
        free( work );
    }

    (void)lowrank;
    (void)Cn;
}

void
core_zlrlr2fr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               pastix_int_t Cm, pastix_int_t Cn,
               pastix_int_t offx, pastix_int_t offy,
               pastix_complex64_t alpha, const pastix_lrblock_t *A,
                                         const pastix_lrblock_t *B,
               pastix_complex64_t beta,        pastix_lrblock_t *C,
               pastix_complex64_t *work, pastix_int_t lwork,
               SolverCblk *fcblk )
{
    pastix_complex64_t *Cptr;
    pastix_int_t        ldcu;
    pastix_lrblock_t    AB;
    pastix_trans_t      trans;
    int                 infomask = 0;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    trans = core_zlrlr2lr( lowrank, transA, transB, M, N, K, A, B, &AB, NULL, -1, &infomask );
    assert( AB.rk != -1 );
    assert( AB.rkmax != -1 );

    if ( AB.rk > 0 ) {
        pastix_int_t ldabv = (trans == PastixNoTrans) ? AB.rkmax : N;

        pastix_cblk_lock( fcblk );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        kernel_trace_start_lvl2( PastixKernelLvl2_LRLR2FR );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                     M, N, AB.rk,
                     CBLAS_SADDR(alpha), AB.u, M,
                                         AB.v, ldabv,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );
        pastix_cblk_unlock( fcblk );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    (void)lowrank;
    (void)Cn;
    (void)work;
    (void)lwork;
}

pastix_trans_t
core_zfrfr2lr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               const pastix_lrblock_t *A,
               const pastix_lrblock_t *B,
               pastix_lrblock_t *AB,
               pastix_complex64_t *work, pastix_int_t lwork,
               int *infomask )
{
    pastix_int_t ldau, ldbu;
    pastix_trans_t transV = PastixNoTrans;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     * Everything is full rank
     */
    if ( K < pastix_imin( M, N ) ) {
        /*
         * Let's build a low-rank matrix of rank K
         */
        AB->rk = K;
        AB->rkmax = K;
        AB->u = A->u;
        AB->v = B->u;
        transV = transB;
    }
    else {
        /*
         * Let's compute the product to form a full-rank matrix of rank
         * pastix_imin( M, N )
         */
        AB->rk = -1;
        AB->rkmax = M;
        if ( lwork < M * N ) {
            work = malloc( M * N * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }
        AB->u = work;
        AB->v = NULL;

        kernel_trace_start_lvl2( PastixKernelLvl2_FRFR2LR );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->u, M );
        kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, K ) );
    }

    (void)lowrank;
    return transV;
}

pastix_trans_t
core_zfrlr2lr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               const pastix_lrblock_t *A,
               const pastix_lrblock_t *B,
               pastix_lrblock_t *AB,
               pastix_complex64_t *work, pastix_int_t lwork,
               int *infomask )
{
    pastix_int_t ldau, ldbu, ldbv;
    pastix_trans_t transV = PastixNoTrans;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( M < B->rk ) {
        /*
         * We are in a similar case to the _Cfr function, and we
         * choose the optimal number of flops.
         */
        pastix_fixdbl_t flops1 = FLOPS_ZGEMM( M, B->rk, K ) + FLOPS_ZGEMM( M, N, B->rk );
        pastix_fixdbl_t flops2 = FLOPS_ZGEMM( K, N, B->rk ) + FLOPS_ZGEMM( M, N, K     );
        pastix_complex64_t *tmp;

        AB->rk    = -1;
        AB->rkmax = M;
        AB->v     = NULL;

        if ( flops1 <= flops2 ) {
            if ( lwork < ( M * B->rk + M * N ) ) {
                work = malloc( (M * B->rk + M * N ) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (A * Bv) * Bu^t
             */
            kernel_trace_start_lvl2( PastixKernelLvl2_FRLR2LR );
            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, B->rk, K,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         M, N, B->rk,
                         CBLAS_SADDR(zone),  tmp,   M,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->u, M );
            kernel_trace_stop_lvl2( flops1 );
        }
        else {
            if ( lwork < ( K * N + M * N ) ) {
                work = malloc( (K * N + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  A * (Bu * Bv^t)^t
             */
            kernel_trace_start_lvl2( PastixKernelLvl2_FRLR2LR );
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         K, N, B->rk,
                         CBLAS_SADDR(zone),  B->u, ldbu,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  K );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             tmp,   K,
                         CBLAS_SADDR(zzero), AB->u, M );

            kernel_trace_stop_lvl2( flops2 );
        }
    }
    else {
        /*
         * B->rk is the smallest rank
         */
        AB->rk    = B->rk;
        AB->rkmax = B->rkmax;
        AB->v     = B->u;
        transV   = transB;

        if ( lwork < ( M * B->rk ) ) {
            work = malloc( M * B->rk * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }
        AB->u = work;

        kernel_trace_start_lvl2( PastixKernelLvl2_FRLR2LR );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->v,  ldbv,
                     CBLAS_SADDR(zzero), AB->u, M );
        kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, B->rk, K ) );
    }

    (void)lowrank;
    return transV;
}


pastix_trans_t
core_zlrfr2lr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               const pastix_lrblock_t *A,
               const pastix_lrblock_t *B,
               pastix_lrblock_t *AB,
               pastix_complex64_t *work, pastix_int_t lwork,
               int *infomask )
{
    pastix_int_t ldau, ldav, ldbu;
    pastix_trans_t transV = PastixNoTrans;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( N < A->rk ) {
        /*
         * We are in a similar case to the _Cfr function, and we
         * choose the optimal number of flops.
         */
        pastix_fixdbl_t flops1 = FLOPS_ZGEMM( A->rk, N, K ) + FLOPS_ZGEMM( M, N, A->rk );
        pastix_fixdbl_t flops2 = FLOPS_ZGEMM( M, K, A->rk ) + FLOPS_ZGEMM( M, N, K     );
        pastix_complex64_t *tmp;

        AB->rk    = -1;
        AB->rkmax = M;
        AB->v     = NULL;

        if ( flops1 <= flops2 ) {
            if ( lwork < (A->rk * N + M * N) ) {
                work = malloc( (A->rk * N + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  Au * (Av^t * B^t)
             */
            kernel_trace_start_lvl2( PastixKernelLvl2_LRFR2LR );
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, K,
                         CBLAS_SADDR(zone),  A->v, ldav,
                                             B->u, ldbu,
                         CBLAS_SADDR(zzero), tmp,  A->rk );

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, N, A->rk,
                         CBLAS_SADDR(zone),  A->u,  ldau,
                                             tmp,   A->rk,
                         CBLAS_SADDR(zzero), AB->u, M );
            kernel_trace_stop_lvl2( flops1 );
        }
        else {
            if ( lwork < (M * K + M * N) ) {
                work = malloc( (M * K + M * N) * sizeof(pastix_complex64_t) );
                *infomask |= PASTIX_LRM3_ALLOCU;
            }

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (Au * Av^t) * B^t
             */
            kernel_trace_start_lvl2( PastixKernelLvl2_LRFR2LR );
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, K, A->rk,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             A->v, ldav,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  tmp,   M,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->u, M );
            kernel_trace_stop_lvl2( flops2 );
        }
    }
    else {
        /*
         * A->rk is the smallest rank
         */
        AB->rk    = A->rk;
        AB->rkmax = A->rk;
        AB->u     = A->u;

        if ( lwork < ( A->rk * N ) ) {
            work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCV;
        }
        AB->v = work;

        kernel_trace_start_lvl2( PastixKernelLvl2_LRFR2LR );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v,  ldav,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->v, AB->rkmax );
        kernel_trace_stop_lvl2( FLOPS_ZGEMM( A->rk, N, K ) );
    }

    (void)lowrank;
    return transV;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_lr_null
 *
 * @brief Compute the product of two low-rank matrices (rank != -1) and returns
 * the result in AB
 *
 *******************************************************************************
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 * @param[in] transA
 *         @arg PastixNoTrans: No transpose, op( A ) = A;
 *         @arg PastixTrans:   Transpose, op( A ) = A';
 *
 * @param[in] transB
 *         @arg PastixNoTrans: No transpose, op( B ) = B;
 *         @arg PastixTrans:   Transpose, op( B ) = B';
 *
 * @param[in] M
 *          The number of rows of the matrix A.
 *
 * @param[in] N
 *          The number of columns of the matrix B.
 *
 * @param[in] K
 *          The number of columns of the matrix A and the number of rows of the
 *          matrix B.
 *
 * @param[in] A
 *          The low-rank representation of the matrix A.
 *
 * @param[in] B
 *          The low-rank representation of the matrix B.
 *
 * @param[out] AB
 *          The low-rank representation of the matrix op(A)op(B).
 *
 *******************************************************************************
 *
 * @return The way the product AB is stored: AB or op(AB).
 *
 *******************************************************************************/
pastix_trans_t
core_zlrlr2lr( const pastix_lr_t *lowrank,
               pastix_trans_t transA, pastix_trans_t transB,
               pastix_int_t M, pastix_int_t N, pastix_int_t K,
               const pastix_lrblock_t *A,
               const pastix_lrblock_t *B,
               pastix_lrblock_t *AB,
               pastix_complex64_t *work, pastix_int_t lwork,
               int *infomask )
{
    pastix_int_t ldau, ldav, ldbu, ldbv;
    int transV = PastixNoTrans;
    pastix_complex64_t *work2;
    pastix_lrblock_t rArB;
    pastix_int_t flops = 0;

    assert( A->rk <= A->rkmax && A->rk > 0 );
    assert( B->rk <= B->rkmax && B->rk > 0 );
    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);

    *infomask = 0;
    ldau = (A->rk == -1) ? A->rkmax : M;
    ldav = A->rkmax;
    ldbu = (B->rk == -1) ? B->rkmax : N;
    ldbv = B->rkmax;

    work2 = malloc( A->rk * B->rk * sizeof(pastix_complex64_t));

    /*
     * Let's compute A * B' = Au Av^h (Bu Bv^h)' with the smallest ws
     */
    kernel_trace_start_lvl2( PastixKernelLvl2_LRLR2LR );
    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                 A->rk, B->rk, K,
                 CBLAS_SADDR(zone),  A->v,  ldav,
                                     B->v,  ldbv,
                 CBLAS_SADDR(zzero), work2, A->rk );
    kernel_trace_stop_lvl2( FLOPS_ZGEMM( A->rk, B->rk, K ) );

    /*
     * Try to compress (Av^h Bv^h')
     */
    lowrank->core_ge2lr( lowrank->tolerance, -1, A->rk, B->rk, work2, A->rk, &rArB );

    kernel_trace_start_lvl2( PastixKernelLvl2_LRLR2LR );

    /*
     * The rank of AB is not smaller than min(rankA, rankB)
     */
    if (rArB.rk == -1){
        if ( A->rk <= B->rk ) {
            /*
             *    ABu = Au
             *    ABv = (Av^h Bv^h') * Bu'
             */
            pastix_complex64_t *work = malloc( A->rk * N * sizeof(pastix_complex64_t));

            //assert( (A->rk * ( N + B->rk )) <= lwork );
            AB->rk    = A->rk;
            AB->rkmax = A->rk;
            AB->u     = A->u;
            AB->v     = work;
            *infomask |= PASTIX_LRM3_ORTHOU;
            *infomask |= PASTIX_LRM3_ALLOCV;

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, B->rk,
                         CBLAS_SADDR(zone),  work2, A->rk,
                                             B->u,  ldbu,
                         CBLAS_SADDR(zzero), AB->v, AB->rkmax );
            flops = FLOPS_ZGEMM( A->rk, N, B->rk );
        }
        else {
            /*
             *    ABu = Au * (Av^h Bv^h')
             *    ABv = Bu'
             */
            pastix_complex64_t *work = malloc( B->rk * M * sizeof(pastix_complex64_t));

            //assert( (B->rk * ( M + A->rk )) <= lwork );
            AB->rk = B->rk;
            AB->rkmax = B->rk;
            AB->u = work;
            AB->v = B->u;
            *infomask |= PASTIX_LRM3_ALLOCU;

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, B->rk, A->rk,
                         CBLAS_SADDR(zone),  A->u,   ldau,
                                             work2,  A->rk,
                         CBLAS_SADDR(zzero), AB->u,  M );
            flops = FLOPS_ZGEMM( M, B->rk, A->rk );

            transV = transB;

            /* free(work); */
        }
    }
    else if (rArB.rk == 0){
        AB->rk    = 0;
        AB->rkmax = 0;
        AB->u = NULL;
        AB->v = NULL;
        *infomask |= PASTIX_LRM3_ORTHOU;
    }
    /**
     * The rank of AB is smaller than min(rankA, rankB)
     */
    else {
        pastix_complex64_t *work = malloc( (M + N) * rArB.rk * sizeof(pastix_complex64_t));

        AB->rk    = rArB.rk;
        AB->rkmax = rArB.rk;
        AB->u = work;
        AB->v = work + M * rArB.rk;
        *infomask |= PASTIX_LRM3_ALLOCU;
        *infomask |= PASTIX_LRM3_ORTHOU;

        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, rArB.rk, A->rk,
                     CBLAS_SADDR(zone),  A->u,   ldau,
                                         rArB.u, A->rk,
                     CBLAS_SADDR(zzero), AB->u,  M );

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     rArB.rk, N, B->rk,
                     CBLAS_SADDR(zone),  rArB.v, rArB.rkmax,
                                         B->u,   ldbu,
                     CBLAS_SADDR(zzero), AB->v,  rArB.rk );

        flops = FLOPS_ZGEMM( M, rArB.rk, A->rk ) + FLOPS_ZGEMM( rArB.rk, N, B->rk );

        /* free(work); */
    }
    core_zlrfree(&rArB);
    free(work2);

    /* Not used for now */
    kernel_trace_stop_lvl2( flops );

    (void)transA;
    (void)work;
    (void)lwork;
    return transV;
}

