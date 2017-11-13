/**
 *
 * @file core_zlrmm_Cnull.c
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
#include "pastix_zlrcores.h"
#include "kernels_trace.h"

void
core_zlrmm_Cnull( core_zlrmm_t *params )
{
    PASTE_CORE_ZLRMM_PARAMS( params );
    pastix_lrblock_t AB;
    pastix_trans_t transV = PastixNoTrans;
    int infomask = 0;
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
                                                core_get_rklimit( Cm, Cn ) ) );
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

    flops = 0.0;
    if ( AB.rk != 0 ) {
        pastix_atomic_lock( lock );
        switch ( C->rk ) {
        case -1:
            /*
             * C became full rank
             */
            flops += core_zlr2fr( params, &AB, transV );
            break;

        case 0:
            /*
             * C is still null
             */
            flops += core_zlr2null( params, &AB, transV, infomask );
            break;

        default:
            /*
             * C is low-rank of rank k
             */
            flops += core_zlr2lr( params, &AB, transV );
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
