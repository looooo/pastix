/**
 *
 * @file gpu_ztrsmsp.c
 *
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX GPU kernel routines
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @date 2023-07-21
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "kernels_trace.h"
#include "pastix_zcuda.h"
#include "pastix_cuda.h"
#include <cublas.h>

static char transstr[3] = { 'N', 'T', 'C' };
static char sidestr[2] = { 'L', 'R' };
static char uplostr[3] = { 'U', 'L', 'A' };
static char diagstr[2] = { 'N', 'U' };

/**
 *******************************************************************************
 *
 * @brief Compute the solve update of a block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cuda_ztrsmsp_2dsub( pastix_side_t          side,
                    pastix_uplo_t          uplo,
                    pastix_trans_t         trans,
                    pastix_diag_t          diag,
                    const SolverCblk      *cblk,
                    pastix_int_t           blok_m,
                    const cuDoubleComplex *A,
                    cuDoubleComplex       *C,
                    cudaStream_t           stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zone  = make_cuDoubleComplex( 1.0, 0.0);
#else
    double zone  =  1.0;
#endif
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset, cblk_m, full_m;
    cuDoubleComplex *Cptr;
    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelTRSMBlok2d );

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    blok   = fblok + blok_m;
    offset = blok->coefind;
    cblk_m = blok->fcblknm;
    full_m = 0;

    cublasSetKernelStream( stream );
    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        Cptr = C + blok->coefind - offset;
        M    = blok_rownbr(blok);
        ldc  = M;

        cublasZtrsm( sidestr[side - PastixLeft],
                     uplostr[uplo - PastixUpper],
                     transstr[trans - PastixNoTrans],
                     diagstr[diag - PastixNonUnit],
                     M, N, zone,
                     A, lda,
                     Cptr, ldc );

        flops += FLOPS_ZTRSM( side, M, N );
        full_m += M;
    }

#if defined(PASTIX_GENERATE_MODEL)
    cudaStreamSynchronize( stream );
#endif
    kernel_trace_stop( blok->inlast, PastixKernelTRSMBlok2d,
                       full_m, N, 0, flops, time );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to one off-diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the A matrix appears on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the A matrix is upper or lower triangular. It has to
 *          be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the A matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the A matrix is unit triangular. It has to be either
 *          PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the correct representation of A.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[inout] C
 *          The pointer to the correct representation of C.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
pastix_fixdbl_t
gpublok_ztrsmsp( pastix_side_t      side,
                 pastix_uplo_t      uplo,
                 pastix_trans_t     trans,
                 pastix_diag_t      diag,
                 const SolverCblk  *cblk,
                 pastix_int_t       blok_m,
                 const void        *A,
                 void              *C,
                 const pastix_lr_t *lowrank,
                 cudaStream_t       stream )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        /* return cuda_ztrsmsp_lrsub( side, uplo, trans, diag, */
        /*                            cblk, blok_m, A, C, lowrank ); */
        assert( 0 );
        return 0;
    }
    else {
        return cuda_ztrsmsp_2dsub( side, uplo, trans, diag,
                                   cblk, blok_m, A, C, stream );
    }

    (void)lowrank;
}
