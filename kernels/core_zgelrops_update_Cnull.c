/**
 *
 * @file core_zgelrops_update_Cnull.c
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

pastix_fixdbl_t
core_zfrfr2null( const pastix_lr_t *lowrank,
                 pastix_trans_t transA, pastix_trans_t transB,
                 pastix_int_t M, pastix_int_t N, pastix_int_t K,
                 pastix_int_t Cm, pastix_int_t Cn,
                 const pastix_lrblock_t *A,
                 const pastix_lrblock_t *B,
                 pastix_lrblock_t *AB,
                 pastix_complex64_t *work, pastix_int_t lwork,
                 int *infomask )
{
    pastix_int_t ldau, ldbu;
    pastix_fixdbl_t flops = 0;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     * Everything is full rank
     */
    if ( ( K < pastix_imin( M, N )        ) &&
         ( K < core_zlr_rklimit( Cm, Cn ) ) )
    {
        /*
         * Let's build a low-rank matrix of rank K
         */
        AB->rk    = K;
        AB->rkmax = K;
        AB->u     = A->u;
        AB->v     = B->u;
        *infomask |= PASTIX_LRM3_TRANSB;
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

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->u, M );
        flops = FLOPS_ZGEMM( M, N, K );
    }

    (void)lowrank;
    return flops;
}

pastix_fixdbl_t
core_zfrlr2null( const pastix_lr_t *lowrank,
                 pastix_trans_t transA, pastix_trans_t transB,
                 pastix_int_t M, pastix_int_t N, pastix_int_t K,
                 pastix_int_t Cm, pastix_int_t Cn,
                 const pastix_lrblock_t *A,
                 const pastix_lrblock_t *B,
                 pastix_lrblock_t *AB,
                 pastix_complex64_t *work, pastix_int_t lwork,
                 int *infomask )
{
    pastix_int_t ldau, ldbu, ldbv;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( ( M < B->rk ) ||
         ( B->rk > core_zlr_rklimit( Cm, Cn ) ) )
    {
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
            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, B->rk, K,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         M, N, B->rk,
                         CBLAS_SADDR(zone),  tmp,  M,
                                             B->u, ldbu,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops1;
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
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         K, N, B->rk,
                         CBLAS_SADDR(zone),  B->u, ldbu,
                                             B->v, ldbv,
                         CBLAS_SADDR(zzero), tmp,  K );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             tmp,  K,
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

        if ( lwork < ( M * B->rk ) ) {
            work = malloc( M * B->rk * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCU;
        }
        AB->u = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u,  ldau,
                                         B->v,  ldbv,
                     CBLAS_SADDR(zzero), AB->u, M );

        flops = FLOPS_ZGEMM( M, B->rk, K );
    }

    (void)lowrank;
    return flops;
}

pastix_fixdbl_t
core_zlrfr2null( const pastix_lr_t *lowrank,
                 pastix_trans_t transA, pastix_trans_t transB,
                 pastix_int_t M, pastix_int_t N, pastix_int_t K,
                 pastix_int_t Cm, pastix_int_t Cn,
                 const pastix_lrblock_t *A,
                 const pastix_lrblock_t *B,
                 pastix_lrblock_t *AB,
                 pastix_complex64_t *work, pastix_int_t lwork,
                 int *infomask )
{
    pastix_int_t ldau, ldav, ldbu;
    pastix_fixdbl_t flops;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( ( N < A->rk ) ||
         ( A->rk > core_zlr_rklimit( Cm, Cn ) ) )
    {
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
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                         A->rk, N, K,
                         CBLAS_SADDR(zone),  A->v, ldav,
                                             B->u, ldbu,
                         CBLAS_SADDR(zzero), tmp,  A->rk );

            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, N, A->rk,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             tmp,  A->rk,
                         CBLAS_SADDR(zzero), AB->u, M );

            flops = flops1;
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
            cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                         M, K, A->rk,
                         CBLAS_SADDR(zone),  A->u, ldau,
                                             A->v, ldav,
                         CBLAS_SADDR(zzero), tmp,  M );

            cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                         M, N, K,
                         CBLAS_SADDR(zone),  tmp,  M,
                                             B->u, ldbu,
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

        if ( lwork < ( A->rk * N ) ) {
            work = malloc( A->rk * N * sizeof(pastix_complex64_t) );
            *infomask |= PASTIX_LRM3_ALLOCV;
        }
        AB->v = work;

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v,  ldav,
                                         B->u,  ldbu,
                     CBLAS_SADDR(zzero), AB->v, AB->rkmax );
        flops = FLOPS_ZGEMM( A->rk, N, K );

        *infomask |= PASTIX_LRM3_ORTHOU;
    }

    (void)lowrank;
    return flops;
}

pastix_fixdbl_t
core_zlrlr2null( const pastix_lr_t *lowrank,
                 pastix_trans_t transA, pastix_trans_t transB,
                 pastix_int_t M, pastix_int_t N, pastix_int_t K,
                 pastix_int_t Cm, pastix_int_t Cn,
                 const pastix_lrblock_t *A,
                 const pastix_lrblock_t *B,
                 pastix_lrblock_t *AB,
                 pastix_complex64_t *work, pastix_int_t lwork,
                 int *infomask )
{
    pastix_fixdbl_t flops;

    flops = core_zlrlr2lr( lowrank, transA, transB, M, N, K, A, B, AB, work, lwork, infomask );

    (void)Cm;
    (void)Cn;
    return flops;
}

void
core_zlrmm_Cnull( const pastix_lr_t *lowrank,
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
    pastix_lrblock_t AB;
    pastix_trans_t transV = PastixNoTrans;
    int infomask = 0;
    double tol = lowrank->tolerance;
    int orthou = -1;
    pastix_fixdbl_t flops;

    assert(transA == PastixNoTrans);
    assert(transB != PastixNoTrans);
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_FRFR2null );
            flops = core_zfrfr2null( lowrank, transA, transB,
                                     M, N, K,
                                     Cm, Cn,
                                     A, B, &AB,
                                     work, lwork, &infomask );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_FRLR2null );
            flops = core_zfrlr2null( lowrank, transA, transB,
                                     M, N, K,
                                     Cm, Cn,
                                     A, B, &AB,
                                     work, lwork, &infomask );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LRFR2null );
            flops = core_zlrfr2null( lowrank, transA, transB,
                                     M, N, K,
                                     Cm, Cn,
                                     A, B, &AB,
                                     work, lwork, &infomask );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LRLR2null );
            flops = core_zlrlr2null( lowrank, transA, transB,
                                     M, N, K,
                                     Cm, Cn,
                                     A, B, &AB,
                                     NULL, -1, &infomask );
            kernel_trace_stop_lvl2( flops );
            assert( AB.rk != -1 );
            assert( AB.rkmax != -1 );
        }
    }

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        transV = transB;
    }

    orthou = infomask & PASTIX_LRM3_ORTHOU;
    assert( orthou != -1 );

    if ( AB.rk != 0 ) {
        pastix_int_t ldabu = M;
        pastix_int_t ldabv = (transV == PastixNoTrans) ? AB.rkmax : N;

        pastix_cblk_lock( fcblk );

        /* C->rk has changed in parallel */
        if ( C->rk == -1 ) {
            pastix_complex64_t *Cfr = C->u;

            /* Add A*B */
            if ( AB.rk == -1 ) {
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );
                core_zgeadd( PastixNoTrans, M, N,
                             alpha, AB.u, M,
                             beta,  Cfr + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( 2. * M * N );
            }
            else {
                kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, ldabu,
                                                 AB.v, ldabv,
                             CBLAS_SADDR(beta),  Cfr + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );
            }
        }
        else if ( C->rk > 0 ) {
            pastix_int_t rmax = core_zlr_rklimit( Cm, Cn );
            pastix_int_t rAB = ( AB.rk == -1 ) ? pastix_imin( M, N ) : AB.rk;

            /*
             * The rank is too big, we need to uncompress/compress C
             */
            if ( (C->rk + rAB) > rmax )
            {
                pastix_complex64_t *Cfr = malloc( Cm * Cn * sizeof(pastix_complex64_t) );

                kernel_trace_start_lvl2( PastixKernelLvl2_LR_UNCOMPRESS );
                cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                             Cm, Cn, C->rk,
                             CBLAS_SADDR(zone),  C->u, Cm,
                                                 C->v, C->rkmax,
                             CBLAS_SADDR(zzero), Cfr,  Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( Cm, Cn, C->rk ) );

                /* Add A*B */
                if ( rAB == -1 ) {
                    kernel_trace_start_lvl2( PastixKernelLvl2_LR_GEMM_PRODUCT );
                    core_zgeadd( PastixNoTrans, M, N,
                                 alpha, AB.u, M,
                                 beta,  Cfr + Cm * offy + offx, Cm );
                    kernel_trace_stop_lvl2( 2. * M * N );
                }
                else {
                    kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
                    cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                 M, N, AB.rk,
                                 CBLAS_SADDR(alpha), AB.u, ldabu,
                                                     AB.v, ldabv,
                                 CBLAS_SADDR(beta),  Cfr + Cm * offy + offx, Cm );
                    kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );
                }

                core_zlrfree(C);

                /* Try to recompress */
                /* TODO: flops+= */
                lowrank->core_ge2lr( tol, -1, Cm, Cn, Cfr, Cm, C );

                free(Cfr);
            }
            else {
                lowrank->core_rradd( lowrank, transV, &alpha,
                                     M,  N,  &AB,
                                     Cm, Cn, C,
                                     offx, offy );
            }
        }
        else {
            /* C is still null rank */
            pastix_int_t rklimit = core_zlr_rklimit( Cm, Cn );

            if ( AB.rk > rklimit ) {
                pastix_complex64_t *work = malloc( Cm * Cn * sizeof(pastix_complex64_t) );

                if ( (M != Cm) || (N != Cn) ) {
                    memset( work, 0, Cm * Cn * sizeof(pastix_complex64_t) );
                }
                kernel_trace_start_lvl2( PastixKernelLvl2_FR_GEMM );
                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, ldabu,
                                                 AB.v, ldabv,
                             CBLAS_SADDR(beta),  work + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );

                lowrank->core_ge2lr( lowrank->tolerance, -1, Cm, Cn, work, Cm, C );

                free(work);
            }
            else {
                if ( !orthou ) {
#if 0
                    if ( AB.rk == -1 ) {
                        pastix_lrblock_t backup;

                        lowrank->core_ge2lr( lowrank->tolerance, rklimit,
                                             M, N, AB.u, ldabu, &backup );

                        core_zlrcpy( lowrank, PastixNoTrans, alpha,
                                     M, N, &backup, Cm, Cn, C,
                                     offx, offy );

                        core_zlrfree( &backup );
                    }
                    else {
                        pastix_complex64_t *ABfr = malloc( M * N * sizeof(pastix_complex64_t) );
                        double norm, eps, normu, normv;

                        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                     M, N, AB.rk,
                                     CBLAS_SADDR(zone), AB.u, ldabu,
                                                        AB.v, ldabv,
                                     CBLAS_SADDR(zzero), ABfr, M );

                        normu = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', M, AB.rk, AB.u, ldabu, NULL );
                        normv = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', AB.rk, N, AB.v, ldabv, NULL );
                        norm  = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', M, N, ABfr, M, NULL );

                        fprintf(stderr, "AB.u= %e, AB.v= %e AB=%e (M= %3ld, N= %3ld, rk= %3d, infomask=%d)\n",
                                normu, normv, norm, M, N, AB.rk, infomask );

                        /* /\* TODO: Check why this fails *\/ */
                        /* if ( transV != PastixNoTrans ) { */
                        /*     fprintf(stderr, "Je suis la !!! \n"); */
                        /* } */
                        /* else { */
                        /*     fprintf(stderr, "Je suis ici !!! \n"); */
                        /* } */
                        //core_zlrorthu( transV, M, N, AB.rk, AB.u, ldabu, AB.v, ldabv );

                        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                     M, N, AB.rk,
                                     CBLAS_SADDR(mzone), AB.u, ldabu,
                                                         AB.v, ldabv,
                                     CBLAS_SADDR(zone),  ABfr, M );

                        normu = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', M, AB.rk, AB.u, ldabu, NULL );
                        normv = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', AB.rk, N, AB.v, ldabv, NULL );
                        norm  = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'M', M, N, ABfr, M, NULL );
                        eps = LAPACKE_dlamch_work( 'e' );

                        fprintf(stderr, "AB.u= %e, AB.v= %e, norm= %e\n", normu, normv, norm );
                        if ( (norm / (eps * AB.rk)) > 10. ) {
                            fprintf(stderr, "norm=%e\n", norm );
                            assert( 0 );
                        }
                        free(ABfr);

                        core_zlrcpy( lowrank, transV, alpha,
                                     M, N, &AB, Cm, Cn, C,
                                     offx, offy );
                    }
#else
                    {
                        pastix_complex64_t *ABfr;
                        pastix_lrblock_t backup;

                        if ( AB.rk > 0 ) {
                            ABfr = malloc( M * N * sizeof(pastix_complex64_t) );

                            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                         M, N, AB.rk,
                                         CBLAS_SADDR(zone), AB.u, ldabu,
                                         AB.v, ldabv,
                                         CBLAS_SADDR(zzero), ABfr, M );
                        }
                        else {
                            ABfr = AB.u;
                        }

                        lowrank->core_ge2lr( lowrank->tolerance, rklimit,
                                             M, N, ABfr, M, &backup );

                        core_zlrcpy( lowrank, PastixNoTrans, alpha,
                                     M, N, &backup, Cm, Cn, C,
                                     offx, offy );

                        core_zlrfree( &backup );
                    }
#endif
                }
                else {
                    core_zlrcpy( lowrank, transV, alpha,
                                 M, N, &AB, Cm, Cn, C,
                                 offx, offy );
                }
            }
        }
        pastix_cblk_unlock( fcblk );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    assert( C->rk <= C->rkmax);
}
