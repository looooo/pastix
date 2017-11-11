/**
 *
 * @file core_zlr2xx.c
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
core_zlr2fr( core_zlrmm_t           *params,
             const pastix_lrblock_t *AB,
             pastix_trans_t          transV )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t ldabu = M;
    pastix_int_t ldabv = (transV == PastixNoTrans) ? AB->rkmax : N;
    pastix_fixdbl_t flops = 0.;
    pastix_complex64_t *Cfr = C->u;
    Cfr += Cm * offy + offx;

    /* TODO: find a suitable name to trace this kind of kernel. */
    if ( AB->rk == -1 ) {
        flops = 2 * M * N;
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_updateCfr );
        core_zgeadd( PastixNoTrans, M, N,
                     alpha, AB->u, ldabu,
                     beta,  Cfr,   Cm );
        kernel_trace_stop_lvl2( flops );
    }
    else {
        flops = FLOPS_ZGEMM( M, N, AB->rk );
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_updateCfr );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                     M, N, AB->rk,
                     CBLAS_SADDR(alpha), AB->u, ldabu,
                                         AB->v, ldabv,
                     CBLAS_SADDR(beta),  Cfr,   Cm );
        kernel_trace_stop_lvl2( flops );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

pastix_fixdbl_t
core_zlr2lr( core_zlrmm_t           *params,
             const pastix_lrblock_t *AB,
             pastix_trans_t          transV )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t rklimit = core_get_rklimit( Cm, Cn );
    pastix_int_t rAB = ( AB->rk == -1 ) ? pastix_imin( M, N ) : AB->rk;
    pastix_int_t ldabu = M;
    pastix_int_t ldabv = (transV == PastixNoTrans) ? AB->rkmax : N;
    pastix_fixdbl_t total_flops = 0.;
    pastix_fixdbl_t flops = 0.;

    /*
     * The rank is too big, we need to uncompress/compress C
     */
    if ( (C->rk + rAB) > rklimit )
    {
        pastix_complex64_t *Cfr, *Coff;
        int allocated = 0;
        if ( (Cfr = core_zlrmm_getws( params, Cm * Cn )) == NULL ) {
            Cfr = malloc( Cm * Cn * sizeof(pastix_complex64_t) );
            allocated = 1;
        }
        Coff = Cfr + Cm * offy + offx;

        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
        cblas_zgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,
                     Cm, Cn, C->rk,
                     CBLAS_SADDR(zone),  C->u, Cm,
                                         C->v, C->rkmax,
                     CBLAS_SADDR(zzero), Cfr,  Cm );
        flops = FLOPS_ZGEMM( Cm, Cn, C->rk );

        /* Add A*B */
        if ( AB->rk == -1 ) {
            core_zgeadd( PastixNoTrans, M, N,
                         alpha, AB->u, M,
                         beta,  Coff,  Cm );
            flops += (2. * M * N);
        }
        else {
            cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                         M, N, AB->rk,
                         CBLAS_SADDR(alpha), AB->u, ldabu,
                                             AB->v, ldabv,
                         CBLAS_SADDR(beta),  Coff,  Cm );
            flops += FLOPS_ZGEMM( M, N, AB->rk );
        }
        kernel_trace_stop_lvl2( flops );
        total_flops += flops;

        /* Try to recompress */
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_recompress );
        core_zlrfree(C); // TODO: Can we give it directly to ge2lr as this
        flops = lowrank->core_ge2lr( lowrank->tolerance, -1, Cm, Cn, Cfr, Cm, C );
        kernel_trace_stop_lvl2_rank( flops, C->rk );
        total_flops += flops;

        if (allocated) {
            free(Cfr);
        }
    }
    /*
     * The rank is not too large, we perform a low-rank update
     */
    else {
        // TODO: Need to get the number of flops
        lowrank->core_rradd( lowrank, transV, &alpha,
                             M,  N,  AB,
                             Cm, Cn, C,
                             offx, offy );
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}

pastix_fixdbl_t
core_zlr2null(core_zlrmm_t           *params,
              const pastix_lrblock_t *AB,
              pastix_trans_t          transV,
              int                     infomask )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_int_t rklimit = core_get_rklimit( Cm, Cn );
    pastix_int_t ldabu = M;
    pastix_int_t ldabv = (transV == PastixNoTrans) ? AB->rkmax : N;
    pastix_fixdbl_t total_flops = 0.;
    pastix_fixdbl_t flops;
    int allocated = 0;

    if ( AB->rk > rklimit ) {
        pastix_complex64_t *Cfr, *Coff;
        if ( (Cfr = core_zlrmm_getws( params, Cm * Cn )) == NULL ) {
            Cfr = malloc( Cm * Cn * sizeof(pastix_complex64_t) );
            allocated = 1;
        }
        Coff = Cfr + Cm * offy + offx;

        /* Set to 0 if contribution smaller than C */
        if ( (M != Cm) || (N != Cn) ) {
            memset( Cfr, 0, Cm * Cn * sizeof(pastix_complex64_t) );
        }

        /* Uncompress the AB product into C */
        flops = FLOPS_ZGEMM( M, N, AB->rk );
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                     M, N, AB->rk,
                     CBLAS_SADDR(alpha), AB->u, ldabu,
                                         AB->v, ldabv,
                     CBLAS_SADDR(beta),  Coff,  Cm );
        kernel_trace_stop_lvl2( flops );
        total_flops += flops;

        /* Try to recompress C */
        kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_recompress );
        flops = lowrank->core_ge2lr( lowrank->tolerance, -1, Cm, Cn, Cfr, Cm, C );
        kernel_trace_stop_lvl2_rank( flops, C->rk );
        total_flops += flops;

        if ( allocated ) {
            free( work );
        }
    }
    else {
        /*
         * Let's chech that AB->u is orthogonal before copying it to C.u
         */
        int orthou = infomask & PASTIX_LRM3_ORTHOU;
        if ( !orthou ) {
            pastix_complex64_t *ABfr;
            pastix_lrblock_t    backup;
            int allocated = 0;

            kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_orthou );

            if ( AB->rk > 0 ) {
                if ( (ABfr = core_zlrmm_getws( params, M * N )) == NULL ) {
                    ABfr = malloc( M * N * sizeof(pastix_complex64_t) );
                    allocated = 1;
                }

                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                             M, N, AB->rk,
                             CBLAS_SADDR(zone),  AB->u, ldabu,
                                                 AB->v, ldabv,
                             CBLAS_SADDR(zzero), ABfr,  M );
                flops = FLOPS_ZGEMM( M, N, AB->rk );
            }
            else {
                ABfr = AB->u;
            }

            flops += lowrank->core_ge2lr( lowrank->tolerance, rklimit,
                                          M, N, ABfr, M, &backup );

            core_zlrcpy( lowrank, PastixNoTrans, alpha,
                         M, N, &backup, Cm, Cn, C,
                         offx, offy );

            kernel_trace_stop_lvl2( flops );
            core_zlrfree( &backup );
            total_flops += flops;

            if ( allocated ) {
                free( ABfr );
            }
        }
        /*
         * AB->u is orthogonal, we directly copy AB->u into C
         */
        else {
            core_zlrcpy( lowrank, transV, alpha,
                         M, N, AB, Cm, Cn, C,
                         offx, offy );
        }
    }

    PASTE_CORE_ZLRMM_VOID;
    return flops;
}
