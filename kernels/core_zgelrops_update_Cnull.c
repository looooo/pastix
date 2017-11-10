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
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t zzero =  0.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

void
core_zlrmm_Cnull( core_zlrmm_t     *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
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
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRFR2null );
            flops = core_zfrfr2lr( params, &AB, &infomask,
                                   pastix_imin( pastix_imin( M, N ),
                                                core_get_rklimit( Cm, Cn )) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_FRLR2null );
            flops = core_zfrlr2lr( params, &AB, &infomask,
                                   pastix_imin( M, core_get_rklimit( Cm, Cn ) ) );
            kernel_trace_stop_lvl2( flops );
        }
    }
    else {
        if ( B->rk == -1 ) {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRFR2null );
            flops = core_zlrfr2lr( params, &AB, &infomask,
                                   pastix_imin( N, core_get_rklimit( Cm, Cn ) ) );
            kernel_trace_stop_lvl2( flops );
        }
        else {
            kernel_trace_start_lvl2( PastixKernelLvl2_LR_LRLR2null );
            flops = core_zlrlr2lr( params, &AB, &infomask );
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

    flops = 0.0;
    if ( AB.rk != 0 ) {
        pastix_int_t ldabu = M;
        pastix_int_t ldabv = (transV == PastixNoTrans) ? AB.rkmax : N;

        pastix_atomic_lock( lock );

        /* C->rk has changed in parallel */
        /* Todo: find a suitable name to trace this kind of kernel.. */
        if ( C->rk == -1 ) {
            pastix_complex64_t *Cfr = C->u;

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
        else if ( C->rk > 0 ) {
            pastix_int_t rmax = core_get_rklimit( Cm, Cn );
            pastix_int_t rAB = ( AB.rk == -1 ) ? pastix_imin( M, N ) : AB.rk;

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
        else {
            /* C is still null rank */
            pastix_int_t rklimit = core_get_rklimit( Cm, Cn );

            if ( AB.rk > rklimit ) {
                pastix_complex64_t *work = malloc( Cm * Cn * sizeof(pastix_complex64_t) );

                if ( (M != Cm) || (N != Cn) ) {
                    memset( work, 0, Cm * Cn * sizeof(pastix_complex64_t) );
                }
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_uncompress );
                cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                             M, N, AB.rk,
                             CBLAS_SADDR(alpha), AB.u, ldabu,
                                                 AB.v, ldabv,
                             CBLAS_SADDR(beta),  work + Cm * offy + offx, Cm );
                kernel_trace_stop_lvl2( FLOPS_ZGEMM( M, N, AB.rk ) );

                kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_recompress );
                flops = lowrank->core_ge2lr( lowrank->tolerance, -1, Cm, Cn, work, Cm, C );
                kernel_trace_stop_lvl2_rank( flops, C->rk );

                free(work);
            }
            else {
                if ( !orthou ) {
                    pastix_complex64_t *ABfr;
                    pastix_lrblock_t    backup;

                    kernel_trace_start_lvl2( PastixKernelLvl2_LR_add2C_orthou );

                    if ( AB.rk > 0 ) {
                        ABfr = malloc( M * N * sizeof(pastix_complex64_t) );

                        cblas_zgemm( CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)transV,
                                     M, N, AB.rk,
                                     CBLAS_SADDR(zone),  AB.u, ldabu,
                                                         AB.v, ldabv,
                                     CBLAS_SADDR(zzero), ABfr, M );
                        flops = FLOPS_ZGEMM( M, N, AB.rk );
                    }
                    else {
                        ABfr = AB.u;
                    }

                    flops += lowrank->core_ge2lr( lowrank->tolerance, rklimit,
                                                  M, N, ABfr, M, &backup );

                    core_zlrcpy( lowrank, PastixNoTrans, alpha,
                                 M, N, &backup, Cm, Cn, C,
                                 offx, offy );

                    kernel_trace_stop_lvl2( flops );
                    core_zlrfree( &backup );
                }
                else {
                    core_zlrcpy( lowrank, transV, alpha,
                                 M, N, &AB, Cm, Cn, C,
                                 offx, offy );
                }
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
