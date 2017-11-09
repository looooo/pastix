/**
 *
 * @file core_zgelrops_update_Cfr.c
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
#include "pastix_zlrcores.h"
#include "kernels_trace.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

static inline pastix_fixdbl_t
core_zfrfr2fr( core_zlrmm_t *params )
{
    pastix_int_t ldau, ldbu, ldcu;
    pastix_complex64_t *Cptr;
    pastix_fixdbl_t flops;
    PASTE_CORE_ZLRMM_PARAMS( params );
    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldcu = Cm;

    Cptr  = C->u;
    Cptr += ldcu * offy + offx;

    pastix_atomic_lock( lock );
    assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */

    /*
     * Everything is full rank we apply directly a GEMM
     */
    cblas_zgemm( CblasColMajor, transA, transB,
                 M, N, K,
                 CBLAS_SADDR(alpha), A->u, ldau,
                                     B->u, ldbu,
                 CBLAS_SADDR(beta),  Cptr, ldcu );
    flops = FLOPS_ZGEMM( M, N, K );

    pastix_atomic_unlock( lock );

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

static inline pastix_fixdbl_t
core_zfrlr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldbu, ldbv, ldcu;
    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( M, B->rk, K ) + FLOPS_ZGEMM( M, N, B->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( K, N, B->rk ) + FLOPS_ZGEMM( M, N, K     );
    pastix_fixdbl_t flops;
    ssize_t wsize;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldbu = (transB == PastixNoTrans) ? K : N;
    ldbv = ( B->rk == -1 ) ? -1 : B->rkmax;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A(M-by-K) * B( N-by-rb x rb-by-K )^t
     */
    if ( flops1 <= flops2 ) {
        wsize = M * B->rk;
        if ( lwork < wsize ) {
            params->work  = work  = realloc( work, wsize * sizeof(pastix_complex64_t) );
            params->lwork = lwork = wsize;
        }

        /*
         *  (A * Bv) * Bu^t
         */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, B->rk, K,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, M );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     M, N, B->rk,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        flops = flops1;
        pastix_atomic_unlock( lock );
    }
    else {
        wsize = K * N;
        if ( lwork < wsize ) {
            params->work  = work  = realloc( work, wsize * sizeof(pastix_complex64_t) );
            params->lwork = lwork = wsize;
        }

        /*
         *  A * (Bu * Bv^t)^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     K, N, B->rk,
                     CBLAS_SADDR(zone),  B->u, ldbu,
                                         B->v, ldbv,
                     CBLAS_SADDR(zzero), work, K );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, K,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops2;
        pastix_atomic_unlock( lock );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

static inline pastix_fixdbl_t
core_zlrfr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t ldau, ldav, ldbu, ldcu;
    pastix_fixdbl_t flops1 = FLOPS_ZGEMM( A->rk, N, K ) + FLOPS_ZGEMM( M, N, A->rk );
    pastix_fixdbl_t flops2 = FLOPS_ZGEMM( M, K, A->rk ) + FLOPS_ZGEMM( M, N, K     );
    pastix_fixdbl_t flops;
    ssize_t wsize;

    ldau = (transA == PastixNoTrans) ? M : K;
    ldav = ( A->rk == -1 ) ? -1 : A->rkmax;
    ldbu = (transB == PastixNoTrans) ? K : N;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    /*
     *  A( M-by-ra x ra-by-K ) * B(N-by-K)^t
     */
    if ( flops1 <= flops2 ) {
        wsize = A->rk * N;
        if ( lwork < wsize ) {
            params->work  = work  = realloc( work, wsize * sizeof(pastix_complex64_t) );
            params->lwork = lwork = wsize;
        }

        /*
         *  Au * (Av^t * B^t)
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transB,
                     A->rk, N, K,
                     CBLAS_SADDR(zone),  A->v, ldav,
                                         B->u, ldbu,
                     CBLAS_SADDR(zzero), work, A->rk );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, N, A->rk,
                     CBLAS_SADDR(alpha), A->u, ldau,
                                         work, A->rk,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops1;
        pastix_atomic_unlock( lock );
    }
    else {
        wsize = M * K;
        if ( lwork < wsize ) {
            params->work  = work  = realloc( work, wsize * sizeof(pastix_complex64_t) );
            params->lwork = lwork = wsize;
        }

        /*
         *  (Au * Av^t) * B^t
         */
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     M, K, A->rk,
                     CBLAS_SADDR(zone),  A->u, ldau,
                                         A->v, ldav,
                     CBLAS_SADDR(zzero), work, M );

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */
        cblas_zgemm( CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                     M, N, K,
                     CBLAS_SADDR(alpha), work, M,
                                         B->u, ldbu,
                     CBLAS_SADDR(beta),  Cptr, ldcu );

        flops = flops2;
        pastix_atomic_unlock( lock );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

static inline pastix_fixdbl_t
core_zlrlr2fr( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_complex64_t *Cptr;
    pastix_int_t        ldcu;
    pastix_lrblock_t    AB;
    pastix_trans_t      trans = PastixNoTrans;
    int                 infomask = 0;
    pastix_fixdbl_t     flops;

    ldcu = Cm;
    Cptr = C->u;
    Cptr += ldcu * offy + offx;

    flops = core_zlrlr2lr( params, &AB, &infomask );
    assert( AB.rk != -1 );
    assert( AB.rkmax != -1 );

    if ( infomask & PASTIX_LRM3_TRANSB ) {
        trans = transB;
    }

    if ( AB.rk > 0 ) {
        pastix_int_t ldabv = (trans == PastixNoTrans) ? AB.rkmax : N;

        pastix_atomic_lock( lock );
        assert( C->rk == -1 ); /* Check that C has not changed due to parallelism */

        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                     M, N, AB.rk,
                     CBLAS_SADDR(alpha), AB.u, M,
                                         AB.v, ldabv,
                     CBLAS_SADDR(beta),  Cptr, ldcu );
        flops = FLOPS_ZGEMM( M, N, AB.rk );
        pastix_atomic_unlock( lock );
    }

    /* Free memory from zlrm3 */
    if ( infomask & PASTIX_LRM3_ALLOCU ) {
        free(AB.u);
    }
    if ( infomask & PASTIX_LRM3_ALLOCV ) {
        free(AB.v);
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

void
core_zlrmm_Cfr( core_zlrmm_t *params )
{
    const pastix_lrblock_t *A = params->A;
    const pastix_lrblock_t *B = params->B;
    pastix_fixdbl_t flops;

    assert( params->transA == PastixNoTrans );
    assert( params->transB != PastixNoTrans );
    assert( A->rk <= A->rkmax && A->rk != 0 );
    assert( B->rk <= B->rkmax && B->rk != 0 );
    assert( params->C->rk == -1 );

    if ( A->rk == -1 ) {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2FR );
            flops = core_zfrfr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2FR );
            flops = core_zfrlr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2FR );
            flops = core_zlrfr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2FR );
            flops = core_zlrlr2fr( params );
            kernel_trace_stop_lvl2( flops );
        }
    }

    assert( params->C->rk == -1 );
}

