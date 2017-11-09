/**
 *
 * @file core_zgelrops_update_Clr.c
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
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

pastix_fixdbl_t
core_zfrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Kmax )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldbu;
    pastix_fixdbl_t flops = 0.0;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     * Everything is full rank
     */
    if ( K < Kmax ) {
        /*
         * Let's build a low-rank matrix of rank K
         */
        AB->rk = K;
        AB->rkmax = K;
        AB->u = A->u;
        AB->v = B->u;
        *infomask |= PASTIX_LRM3_TRANSB;
    }
    else {
        /*
         * Let's compute the product to form a full-rank matrix of rank
         * pastix_imin( M, N )
         */
        work = core_zlrmm_resize_ws( params, M * N );

        AB->rk = -1;
        AB->rkmax = M;
        AB->u = work;
        AB->v = NULL;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->u, M );
        flops = FLOPS_ZGEMM( M, N, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

pastix_fixdbl_t
core_zfrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Brkmin )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldbu, ldbv;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( B->rk > Brkmin ) {
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
            work = core_zlrmm_resize_ws( params, M * B->rk + M * N );

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (A * Bv) * Bu^t
             */
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

            flops = flops1;
        }
        else {
            work = core_zlrmm_resize_ws( params, K * N + M * N );

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  A * (Bu * Bv^t)^t
             */
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

            flops = flops2;
        }
    }
    else {
        /*
         * B->rk is the smallest rank
         */
        AB->rk    = B->rk;
        AB->rkmax = B->rkmax;
        AB->v     = B->u;
        *infomask |= PASTIX_LRM3_TRANSB;

        work = core_zlrmm_resize_ws( params, M * B->rk );
        AB->u = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->v,  ldbv,
                     CBLAS_SADDR(zzero), AB->u, M );
        flops = FLOPS_ZGEMM( M, B->rk, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

pastix_fixdbl_t
core_zlrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Arkmin )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldav, ldbu;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( A->rk > Arkmin ) {
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
            work = core_zlrmm_resize_ws( params, A->rk * N + M * N );

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  Au * (Av^t * B^t)
             */
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

            flops = flops1;
        }
        else {
            work = core_zlrmm_resize_ws( params, M * K + M * N );

            /* AB->u will be destroyed later */
            AB->u = work;
            tmp   = work + M * N;

            /*
             *  (Au * Av^t) * B^t
             */
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

            flops = flops2;
        }
    }
    else {
        /*
         * A->rk is the smallest rank
         */
        AB->rk    = A->rk;
        AB->rkmax = A->rk;
        AB->u     = A->u;

        work = core_zlrmm_resize_ws( params, A->rk * N );
        AB->v = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v,  ldav,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->v, AB->rkmax );

        flops = FLOPS_ZGEMM( A->rk, N, K );
    }

    PASTE_CORE_ZLRMM_VOID;
    (void)infomask;
    return flops;
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
pastix_fixdbl_t
core_zlrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldau, ldav, ldbu, ldbv;
    pastix_complex64_t *work2;
    pastix_lrblock_t rArB;
    pastix_fixdbl_t flops = 0.0;

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
    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                 A->rk, B->rk, K,
                 CBLAS_SADDR(zone),  A->v,  ldav,
                                     B->v,  ldbv,
                 CBLAS_SADDR(zzero), work2, A->rk );
    flops = FLOPS_ZGEMM( A->rk, B->rk, K );

    /*
     * Try to compress (Av^h Bv^h')
     */

    flops += lowrank->core_ge2lr( lowrank->tolerance, -1, A->rk, B->rk, work2, A->rk, &rArB );

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
            flops += FLOPS_ZGEMM( A->rk, N, B->rk );
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
            flops += FLOPS_ZGEMM( M, B->rk, A->rk );

            *infomask |= PASTIX_LRM3_TRANSB;
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

        flops += FLOPS_ZGEMM( M, rArB.rk, A->rk ) + FLOPS_ZGEMM( rArB.rk, N, B->rk );

        /* free(work); */
    }
    core_zlrfree(&rArB);
    free(work2);

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

void
core_zlrmm_Clr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_lrblock_t AB;
    pastix_trans_t transV = PastixNoTrans;
    int infomask = 0;
    double tol = lowrank->tolerance;
    pastix_fixdbl_t flops;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2LR );
            flops = core_zfrfr2lr( params, &AB, &infomask,
                                   pastix_imin( M, N ) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2LR );
            flops = core_zfrlr2lr( params, &AB, &infomask, M );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2LR );
            flops = core_zlrfr2lr( params, &AB, &infomask, N );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2LR );
            flops = core_zlrlr2lr( params, &AB, &infomask );
            kernel_trace_stop_lvl2( flops );

            assert( AB.rk != -1 );
            assert( AB.rkmax != -1 );
        }
    }

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        transV = transB;
    }

    if ( AB.rk != 0 ) {
        pastix_atomic_lock( lock );

        /* C->rk has changed in parallel */
        /* Todo: find a suitable name to trace this kind of kernel.. */
        if ( C->rk == -1 ) {
            pastix_complex64_t *Cfr = C->u;
            pastix_int_t ldabu = M;
            pastix_int_t ldabv = (transV == PastixNoTrans) ? AB.rkmax : N;

            /* Add A*B */
            if ( AB.rk == -1 ) {
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_updateCfr );
                core_zgeadd( PastixNoTrans, M, N,
                             alpha, AB.u, M,
                             beta,  Cfr + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( 2. * M * N );
            }
            else {
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_updateCfr );
                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, ldabu,
                             AB.v, ldabv,
                             CBLAS_SADDR(beta),  Cfr + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );
            }
        }
        else {
            pastix_int_t rmax = core_zlr_rklimit( Cm, Cn );
            pastix_int_t rAB = ( AB.rk == -1 ) ? pastix_imin( M, N ) : AB.rk;
            pastix_int_t ldabu = M;
            pastix_int_t ldabv = (transV == PastixNoTrans) ? AB.rkmax : N;

            /*
             * The rank is too big, we need to uncompress/compress C
             */
            if ( (C->rk + rAB) > rmax )
            {
                pastix_complex64_t *Cfr = malloc( Cm * Cn * sizeof(pastix_complex64_t) );

                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                             Cm, Cn, C->rk,
                             CBLAS_SADDR(zone),  C->u, Cm,
                                                 C->v, C->rkmax,
                             CBLAS_SADDR(zzero), Cfr,  Cm );
                flops = FLOPS_ZGEMM( Cm, Cn, C->rk );

                /* Add A*B */
                if ( AB.rk == -1 ) {
                    core_zgeadd( PastixNoTrans, M, N,
                                 alpha, AB.u, M,
                                 beta,  Cfr + Cm * offy + offx, Cm );
                    flops += (2. * M * N);
                }
                else {
                    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                 M, N, AB.rk,
                                 CBLAS_SADDR(alpha), AB.u, ldabu,
                                                     AB.v, ldabv,
                                 CBLAS_SADDR(beta),  Cfr + Cm * offy + offx, Cm );
                    flops += FLOPS_ZGEMM( M, N, AB.rk );
                }
                kernel_trace_stop_lvl2( flops );

                core_zlrfree(C);

                /* Try to recompress */
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_recompress );
                flops = lowrank->core_ge2lr( tol, -1, Cm, Cn, Cfr, Cm, C );
                kernel_trace_stop_lvl2_rank( flops, C->rk );

                free(Cfr);
            }
            else {
                lowrank->core_rradd( lowrank, transV, &alpha,
                                     M,  N,  &AB,
                                     Cm, Cn, C,
                                     offx, offy );
            }
        }
        pastix_atomic_unlock( lock );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    assert( C->rk <= C->rkmax);
    PASTE_CORE_ZLRMM_VOID;
}
