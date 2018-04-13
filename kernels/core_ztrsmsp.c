/**
 *
 * @file core_ztrsmsp.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static pastix_complex64_t  zone =  1.0;
static pastix_complex64_t zzero =  0.0;
static pastix_complex64_t mzone = -1.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Apply all the trsm updates on a panel stored in 1D layout.
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
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
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
 *******************************************************************************/
static inline void
core_ztrsmsp_1d( pastix_side_t             side,
                 pastix_uplo_t             uplo,
                 pastix_trans_t            trans,
                 pastix_diag_t             diag,
                 SolverCblk               *cblk,
                 const pastix_complex64_t *A,
                 pastix_complex64_t       *C )
{
    SolverBlok *fblok;
    pastix_int_t M, N, lda;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    lda   = cblk->stride;
    fblok = cblk->fblokptr;  /* The diagonal block */

    /* vertical dimension */
    M = lda - N;

    /* if there is an extra-diagonal bloc in column block */
    assert( fblok + 1 < cblk[1].fblokptr );
    assert( blok_rownbr( fblok) == N );
    assert(!(cblk->cblktype & CBLK_LAYOUT_2D));

    /* first extra-diagonal bloc in column block address */
    C = C + fblok[1].coefind;

    kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
    cblas_ztrsm(CblasColMajor,
                (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                M, N,
                CBLAS_SADDR(zone), A, lda,
                                   C, lda);
    kernel_trace_stop_lvl2( FLOPS_ZTRSM( side, M, N ) );
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in 2D.
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
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
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
 *******************************************************************************/
static inline void
core_ztrsmsp_2d( pastix_side_t             side,
                 pastix_uplo_t             uplo,
                 pastix_trans_t            trans,
                 pastix_diag_t             diag,
                 SolverCblk               *cblk,
                 const pastix_complex64_t *A,
                 pastix_complex64_t       *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc;
    pastix_complex64_t *blokC;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    for (blok=fblok+1; blok<lblok; blok++) {

        blokC = C + blok->coefind;
        M   = blok_rownbr(blok);
        ldc = M;

        kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
        cblas_ztrsm(CblasColMajor,
                    (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       blokC, ldc);
        kernel_trace_stop_lvl2( FLOPS_ZTRSM( side, M, N ) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Computes the updates associated to one off-diagonal block between two
 * cblk stored in low-rank format.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************
 *
 * @return  The number of flops performed
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
core_ztrsmsp_lr( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                 SolverCblk *cblk, const pastix_lr_t *lowrank )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda;
    pastix_lrblock_t *lrA, *lrC;
    pastix_complex64_t *A;

    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t flops_c;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( lrA->rk == -1 );
    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_COMPRESSED );
    assert( cblk->cblktype & CBLK_LAYOUT_2D  );

    for (blok=fblok+1; blok<lblok; blok++) {

        lrC = blok->LRblock + coef;

        /* Try to compress the block: compress_end version */
        if ( lowrank->compress_when == PastixCompressWhenEnd )
        {

            M = blok_rownbr(blok);
            if ( ( N > lowrank->compress_min_width  ) &&
                 ( M > lowrank->compress_min_height ) )
            {
                pastix_lrblock_t C;

                kernel_trace_start_lvl2( PastixKernelLvl2_LR_init_compress );
                flops_c = lowrank->core_ge2lr( lowrank->tolerance, -1,
                                             M, N, lrC->u, M, &C );
                kernel_trace_stop_lvl2_rank( flops_c, C.rk );
                flops += flops_c;

                core_zlrfree(lrC);
                memcpy( lrC, &C, sizeof(pastix_lrblock_t) );
            }
        }

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                kernel_trace_start_lvl2( PastixKernelLvl2_LR_TRSM );
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);
                flops_c = FLOPS_ZTRSM( side, lrC->rk, N );
                kernel_trace_stop_lvl2( flops_c );
                flops += flops_c;
            }
            else {
                M = blok_rownbr(blok);
                kernel_trace_start_lvl2( PastixKernelLvl2_FR_TRSM );
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);
                flops_c = FLOPS_ZTRSM( side, M, N );
                kernel_trace_stop_lvl2( flops_c );
                flops += flops_c;
            }
        }
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to a column of off-diagonal blocks.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          Specify whether we work with the lower matrix, or the upper matrix.
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
void
cpucblk_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C,
                 const pastix_lr_t        *lowrank )
{
    if (  cblk[0].fblokptr + 1 < cblk[1].fblokptr )
    {
        pastix_ktype_t ktype = -1;
        pastix_fixdbl_t time, flops = 0.0;
        pastix_int_t n = cblk_colnbr( cblk );
        pastix_int_t m = cblk->stride - n;

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            ktype = PastixKernelTRSMCblkLR;
            time  = kernel_trace_start( ktype );

            flops = core_ztrsmsp_lr( coef, side, uplo, trans, diag,
                                     cblk, lowrank );
        }
        else {
            if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
                ktype = PastixKernelTRSMCblk2d;
                time  = kernel_trace_start( ktype );

                core_ztrsmsp_2d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
            else {
                ktype = PastixKernelTRSMCblk1d;
                time  = kernel_trace_start( ktype );

                core_ztrsmsp_1d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
            flops = FLOPS_ZTRSM( PastixRight, m, n );
        }

        kernel_trace_stop( ktype, m, n, 0, flops, time );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in 2D.
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
 *          The cblk structure to which block belongs to. The A and C pointers
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
 *******************************************************************************/
static inline void
core_ztrsmsp_2dsub( pastix_side_t side, pastix_uplo_t uplo,
                    pastix_trans_t trans, pastix_diag_t diag,
                          SolverCblk         *cblk,
                          pastix_int_t        blok_m,
                    const pastix_complex64_t *A,
                          pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset, cblk_m, full_m;
    pastix_complex64_t *Cptr;
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

    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        Cptr = C + blok->coefind - offset;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm( CblasColMajor,
                     (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                     M, N,
                     CBLAS_SADDR(zone), A, lda,
                                        Cptr, ldc );

        flops += FLOPS_ZTRSM( side, M, N );
        full_m += M;
    }

    kernel_trace_stop( PastixKernelTRSMBlok2d,
                       full_m, N, 0, flops, time );
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact_null
 *
 * @brief Compute the updates associated to one off-diagonal block between two
 * cblk stored in low-rank format.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right in the
 *          equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************/
static inline void
core_ztrsmsp_lrsub( pastix_coefside_t   coef,
                    pastix_side_t       side,
                    pastix_uplo_t       uplo,
                    pastix_trans_t      trans,
                    pastix_diag_t       diag,
                    SolverCblk         *cblk,
                    pastix_int_t        blok_m,
                    const pastix_lr_t  *lowrank )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, cblk_m, full_m, full_n;
    pastix_complex64_t *A;
    pastix_lrblock_t *lrA, *lrC;
    pastix_fixdbl_t flops = 0.0;
    pastix_fixdbl_t time = kernel_trace_start( PastixKernelTRSMBlokLR );

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( cblk->cblktype & CBLK_COMPRESSED );
    assert( cblk->cblktype & CBLK_LAYOUT_2D  );

    assert( blok_rownbr(fblok) == N );
    assert( lrA->rk == -1 );

    blok   = fblok + blok_m;
    cblk_m = blok->fcblknm;
    full_m = 0;
    full_n = 0;

    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        M = blok_rownbr(blok);
        lrC = blok->LRblock + coef;

        /* Try to compress the block: compress_end version */
        if ( lowrank->compress_when == PastixCompressWhenEnd )
        {
            if ( ( N > lowrank->compress_min_width ) &&
                 ( M > lowrank->compress_min_height ) )
            {
                pastix_fixdbl_t  flops;
                pastix_lrblock_t C;

                kernel_trace_start_lvl2( PastixKernelLvl2_LR_init_compress );
                flops = lowrank->core_ge2lr( lowrank->tolerance, -1,
                                             M, N, lrC->u, M, &C );
                kernel_trace_stop_lvl2_rank( flops, C.rk );

                core_zlrfree(lrC);
                lrC->u = C.u;
                lrC->v = C.v;
                lrC->rk = C.rk;
                lrC->rkmax = C.rkmax;
            }
        }

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);

                flops += FLOPS_ZTRSM( side, lrC->rk, N );
                full_n += lrC->rk;
            }
            else {
                cblas_ztrsm(CblasColMajor,
                            (CBLAS_SIDE)side, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);

                flops += FLOPS_ZTRSM( side, M, N );
                full_n += M;
            }
        }
        full_m += M;
    }

    kernel_trace_stop( PastixKernelTRSMBlokLR,
                       full_m, N, full_n, flops, time );
}

/**
 *******************************************************************************
 *
 * @brief Compute the updates associated to one off-diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          - PastixLCoef, use the lower part of the off-diagonal blocks.
 *          - PastixUCoef, use the upper part of the off-diagonal blocks
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
void
cpublok_ztrsmsp( pastix_coefside_t coef, pastix_side_t side, pastix_uplo_t uplo,
                 pastix_trans_t trans, pastix_diag_t diag,
                       SolverCblk         *cblk,
                       pastix_int_t        blok_m,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C,
                 const pastix_lr_t        *lowrank )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        core_ztrsmsp_lrsub( coef, side, uplo, trans, diag,
                            cblk, blok_m, lowrank );
    }
    else {
        core_ztrsmsp_2dsub( side, uplo, trans, diag,
                            cblk, blok_m, A, C );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply the solve related to one cblk to all the right hand side.
 *
 *******************************************************************************
 *
 * @param[in] mode
 *          Specify whether the schur complement and interface are applied to
 *          the right-hand-side. It has to be either PastixSolvModeLocal,
 *          PastixSolvModeInterface or PastixSolvModeSchur.
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] datacode
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] nrhs
 *          The number of right hand side.
 *
 * @param[inout] b
 *          The pointer to vectors of the right hand side
 *
 * @param[in] ldb
 *          The leading dimension of b
 *
 *******************************************************************************/
void
solve_ztrsmsp( pastix_solv_mode_t  mode,
               pastix_side_t       side,
               pastix_uplo_t       uplo,
               pastix_trans_t      trans,
               pastix_diag_t       diag,
               SolverMatrix       *datacode,
               SolverCblk         *cblk,
               int                 nrhs,
               pastix_complex64_t *b,
               int                 ldb )
{
    SolverCblk *fcbk;
    SolverBlok *blok;
    pastix_complex64_t *A, *tmp;
    pastix_int_t j, tempm, tempn, lda;
    pastix_lrblock_t *lrA;

    tempn = cblk->lcolnum - cblk->fcolnum + 1;
    lda = (cblk->cblktype & CBLK_LAYOUT_2D) ? blok_rownbr( cblk->fblokptr ) : cblk->stride;

    /*
     *  Left / Upper / NoTrans (Backward)
     */
    if (side == PastixLeft) {
        if (uplo == PastixUpper) {
            /*
             * If cblk is in the schur complement, all brow blocks are in
             * the interface.  Thus, it doesn't generate any update in local
             * mode, and we know that we are at least in interface mode
             * after this test.
             */
            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
                for (j = cblk[0].brownum; j < cblk[1].brownum; j++ ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;

                    if ( fcbk->cblktype & CBLK_IN_SCHUR ) {
                        break;
                    }
                    pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                }
                return;
            }

            if ( cblk->cblktype & CBLK_COMPRESSED ) {
                A = (pastix_complex64_t*)(cblk->fblokptr->LRblock[1].u);
                assert( cblk->fblokptr->LRblock[1].rkmax == lda );
            }
            else {
                A = (pastix_complex64_t*)(cblk->ucoeftab);
            }

            /*  We store U^t, so we swap uplo and trans */
            if (trans == PastixNoTrans) {

                if ( !(cblk->cblktype & CBLK_IN_SCHUR) || (mode == PastixSolvModeSchur) ) {
                    /* Solve the diagonal block */
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasTrans, (CBLAS_DIAG)diag,
                        tempn, nrhs,
                        CBLAS_SADDR(zone), A, lda,
                                           b + cblk->lcolidx, ldb );
                }

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;

                    if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeInterface) ) {
                        continue;
                    }

                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->ucoeftab);
                    lda = (fcbk->cblktype & CBLK_LAYOUT_2D) ? tempn : fcbk->stride;

                    pastix_cblk_lock( fcbk );
                    if ( fcbk->cblktype & CBLK_COMPRESSED ) {
                        lrA = blok->LRblock + 1;

                        switch (lrA->rk){
                        case 0:
                            break;
                        case -1:
                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                tempm, nrhs, tempn,
                                CBLAS_SADDR(mzone), lrA->u, tempn,
                                                    b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                                break;
                        default:
                            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);
                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                lrA->rk, nrhs, tempn,
                                CBLAS_SADDR(zone),  lrA->u, tempn,
                                                    b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zzero), tmp,  lrA->rk );

                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                tempm, nrhs, lrA->rk,
                                CBLAS_SADDR(mzone), lrA->v, lrA->rkmax,
                                                    tmp, lrA->rk,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                            memFree_null(tmp);
                            break;
                        }
                    }
                    else{
                        cblas_zgemm(
                            CblasColMajor, CblasTrans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                                b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                    }
                    pastix_cblk_unlock( fcbk );
                    pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                }
            }
            /*
             *  Left / Upper / [Conj]Trans (Forward)
             */
            else {
                if ( cblk->cblktype & CBLK_IN_SCHUR )
                    return;

                assert(0 /* Not implemented */);
            }
        }
        else {
            A = (pastix_complex64_t*)(cblk->lcoeftab);

            /*
             *  Left / Lower / NoTrans (Forward)
             */
            if (trans == PastixNoTrans) {

                if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) )
                    return;

                /* In sequential */
                assert( cblk->fcolnum == cblk->lcolidx );

                if ( cblk->cblktype & CBLK_COMPRESSED ) {

                    /* Solve the diagonal block */
                    lrA = cblk->fblokptr->LRblock;
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasNoTrans, (CBLAS_DIAG)diag,
                        tempn, nrhs, CBLAS_SADDR(zone),
                        lrA->u, tempn,
                        b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
                        fcbk  = datacode->cblktab + blok->fcblknm;

                        if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
                            return;
                        }

                        tempm = blok->lrownum - blok->frownum + 1;
                        lrA   = blok->LRblock;

                        assert( blok->frownum >= fcbk->fcolnum );
                        assert( tempm <= (fcbk->lcolnum - fcbk->fcolnum + 1));

                        pastix_cblk_lock( fcbk );
                        switch (lrA->rk){
                        case 0:
                            break;
                        case -1:
                            assert( lrA->rkmax == tempm );
                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                tempm, nrhs, tempn,
                                CBLAS_SADDR(mzone), lrA->u, tempm,
                                                    b + cblk->lcolidx, ldb,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );
                            break;
                        default:
                            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);

                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                lrA->rk, nrhs, tempn,
                                CBLAS_SADDR(zone), lrA->v, lrA->rkmax,
                                                   b + cblk->lcolidx, ldb,
                                CBLAS_SADDR(zzero), tmp, lrA->rk);

                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                tempm, nrhs, lrA->rk,
                                CBLAS_SADDR(mzone), lrA->u, tempm,
                                                    tmp, lrA->rk,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );

                            memFree_null(tmp);
                        }
                        pastix_cblk_unlock( fcbk );
                        pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                    }
                }
                else {
                    /* Solve the diagonal block */
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasNoTrans, (CBLAS_DIAG)diag,
                        tempn, nrhs,
                        CBLAS_SADDR(zone), A, lda,
                                           b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
                        fcbk  = datacode->cblktab + blok->fcblknm;

                        if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
                            return;
                        }

                        tempm = blok->lrownum - blok->frownum + 1;

                        assert( blok->frownum >= fcbk->fcolnum );
                        assert( tempm <= (fcbk->lcolnum - fcbk->fcolnum + 1));

                        lda = (cblk->cblktype & CBLK_LAYOUT_2D) ? tempm : cblk->stride;

                        pastix_cblk_lock( fcbk );
                        cblas_zgemm(
                            CblasColMajor, CblasNoTrans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                                b + cblk->lcolidx, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );

                        pastix_cblk_unlock( fcbk );
                        pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                    }

                }
            }
            /*
             *  Left / Lower / [Conj]Trans (Backward)
             */
            else {
                /*
                 * If cblk is in the schur complement, all brow blocks are in
                 * the interface.  Thus, it doesn't generate any update in local
                 * mode, and we know that we are at least in interface mode
                 * after this test.
                 */
                if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
                    for (j = cblk[0].brownum; j < cblk[1].brownum; j++ ) {
                        blok = datacode->bloktab + datacode->browtab[j];
                        fcbk = datacode->cblktab + blok->lcblknm;

                        if ( fcbk->cblktype & CBLK_IN_SCHUR ) {
                            break;
                        }
                        pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                    }
                    return;
                }

                /* For cblk in the schur, only apply the update on the non Schur part */
                if ( cblk->cblktype & CBLK_COMPRESSED ) {
                    A = (pastix_complex64_t*)(cblk->fblokptr->LRblock[0].u);
                    assert( cblk->fblokptr->LRblock[0].rkmax == lda );
                } else {
                    A = (pastix_complex64_t*)(cblk->lcoeftab);
                }

                /* Solve the diagonal block */
                if ( !(cblk->cblktype & CBLK_IN_SCHUR) || (mode == PastixSolvModeSchur) ) {
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
                        tempn, nrhs,
                        CBLAS_SADDR(zone), A, lda,
                                           b + cblk->lcolidx, ldb );
                }

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;
                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->lcoeftab);
                    lda = (fcbk->cblktype & CBLK_LAYOUT_2D) ? tempn : fcbk->stride;

                    if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeInterface) ) {
                        continue;
                    }

                    pastix_cblk_lock( fcbk );
                    if ( fcbk->cblktype & CBLK_COMPRESSED ) {
                        lrA = blok->LRblock;

                        switch (lrA->rk){
                        case 0:
                            break;
                        case -1:
                            cblas_zgemm(
                                CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                tempm, nrhs, tempn,
                                CBLAS_SADDR(mzone), lrA->u, tempn,
                                b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                            break;
                        default:
                            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);
                            /*
                             * TODO : check if trans (if conj is already
                             * applied) or [conj]trans based on case as for v.
                             */
                            cblas_zgemm(
                                CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                lrA->rk, nrhs, tempn,
                                CBLAS_SADDR(zone),  lrA->u, tempn,
                                                    b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zzero), tmp,  lrA->rk );

                            cblas_zgemm(
                                CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                                tempm, nrhs, lrA->rk,
                                CBLAS_SADDR(mzone), lrA->v, lrA->rkmax,
                                                    tmp, lrA->rk,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                            memFree_null(tmp);
                            break;
                        }
                    }
                    else{
                        cblas_zgemm(
                            CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                                b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                    }
                    pastix_cblk_unlock( fcbk );
                    pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
                }
            }
        }
    }
    /**
     * Right
     */
    else {
        assert(0 /* Not implemented */);
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply the diagonal solve related to one cblk to all the right hand side.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          The cblk structure to which diagonal block belongs to.
 *
 * @param[in] nrhs
 *          The number of right hand side
 *
 * @param[inout] b
 *          The pointer to vectors of the right hand side
 *
 * @param[in] ldb
 *          The leading dimension of b
 *
 * @param[inout] work
 *          Workspace to temporarily store the diagonal when multiple RHS are
 *          involved. Might be set to NULL for internal allocation on need.
 *
 *******************************************************************************/
void
solve_zdiag( SolverCblk         *cblk,
             int                 nrhs,
             pastix_complex64_t *b,
             int                 ldb,
             pastix_complex64_t *work )
{
    pastix_complex64_t *A, *tmp;
    pastix_int_t k, j, tempn, lda;

    tempn = cblk->lcolnum - cblk->fcolnum + 1;
    lda = (cblk->cblktype & CBLK_LAYOUT_2D) ? tempn : cblk->stride;
    assert( blok_rownbr( cblk->fblokptr ) == tempn );

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        A = (pastix_complex64_t*)(cblk->fblokptr->LRblock[0].u);
        assert( cblk->fblokptr->LRblock[0].rkmax == lda );
    }
    else {
        A = (pastix_complex64_t*)(cblk->lcoeftab);
    }

    /* Add shift for diagonal elements */
    lda++;

    if( nrhs == 1 ) {
        for (j=0; j<tempn; j++, b++, A+=lda) {
            *b = (*b) / (*A);
        }
    }
    else {
        /* Copy the diagonal to a temporary buffer */
        tmp = work;
        if ( work == NULL ) {
            MALLOC_INTERN( tmp, tempn, pastix_complex64_t );
        }
        cblas_zcopy( tempn, A, lda, tmp, 1 );

        /* Compute */
        for (k=0; k<nrhs; k++, b+=ldb)
        {
            for (j=0; j<tempn; j++) {
                b[j] /= tmp[j];
            }
        }

        if ( work == NULL ) {
            memFree_null(tmp);
        }
    }
}
