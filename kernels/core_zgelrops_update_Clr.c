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

void
core_zlrmm_Clr( core_zlrmm_t *params )
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
