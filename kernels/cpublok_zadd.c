/**
 *
 * @file cpublok_zadd.c
 *
 * Precision dependent routines to add different bloks.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-12-01
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "kernels_trace.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 * @param[in] blokA_m
 *          Index of the first off-diagonal block in cblk, that is used for A.
 *
 * @param[in] blokB_m
 *          Index of the first off-diagonal block in cblk, that is used for B.
 *
 * @param[inout] A
 *          The pointer to the coeftab of the blok.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          blok.ucoeftab otherwise. Must be of size blok.stride -by- blok.width
 *
 * @param[in] lrB
 *          Pointer to the low-rank representation of the column block B.
 *          Must be followed by the low-rank representation of the following blocks.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpublok_zadd_frlr( pastix_int_t              alpha,
                   const SolverCblk         *cblkA,
                   SolverCblk               *cblkB,
                   pastix_int_t              blokA_m,
                   pastix_int_t              blokB_m,
                   const pastix_complex64_t *A,
                   pastix_lrblock_t         *lrB,
                   pastix_complex64_t       *work,
                   pastix_int_t              lwork,
                   const pastix_lr_t        *lowrank )
{
    /*
     * We can start both blocks at the m^th in the cblk, as the cblkB always
     * holds at least as many blocks as A
     */
    const SolverBlok   *blokA  = cblkA->fblokptr + blokA_m;
    const SolverBlok   *blokB  = cblkB->fblokptr + blokB_m;
    const SolverBlok   *lblokA = cblkA[1].fblokptr;
    const SolverBlok   *lblokB = cblkB[1].fblokptr;
    pastix_fixdbl_t     flops = 0.;
    core_zlrmm_t        params;
    pastix_lrblock_t    lrA;
    size_t              offsetA;

    assert( !(cblkA->cblktype & CBLK_COMPRESSED) );
    assert(   cblkB->cblktype & CBLK_COMPRESSED  );
    assert(   cblkA->cblktype & CBLK_LAYOUT_2D   );

    assert( A != NULL );

    params.lowrank = lowrank;
    params.transA  = PastixNoTrans; /* Unused */
    params.transB  = PastixNoTrans; /* Unused */
    params.K       = -1;            /* Unused */
    params.alpha   = alpha;
    params.A       = NULL;          /* Unused */
    params.B       = NULL;          /* Unused */
    params.beta    = 1.0;
    params.work    = work;
    params.lwork   = lwork;
    params.lwused  = 0;
    params.lock    = &(cblkB->lock);

    /* Dimensions on N */
    params.N    = cblk_colnbr( cblkA );
    params.Cn   = cblk_colnbr( cblkB );
    params.offy = cblkA->fcolnum - cblkB->fcolnum;

    lrA.rk = -1;
    lrA.v  = NULL;

    offsetA = blokA->coefind;
    do {
        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++; lrB++;
        }
        assert( blokA->fcblknm == blokB->fcblknm );
        assert( is_block_inside_fblock( blokA, blokB ) );
        assert( blokB <= lblokB );

        lrA.u     = (pastix_complex64_t*)A + blokA->coefind - offsetA;
        lrA.rkmax = blok_rownbr( blokA );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = lrB;

        flops += core_zlradd( &params, &lrA,
                              PastixNoTrans, 0 );
        blokA++;
    }
    while ( ( blokA < lblokA ) &&
            ( blokA[-1].fcblknm == blokA[0].fcblknm ) &&
            ( blokA[-1].lcblknm == blokA[0].lcblknm ) );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 * @param[in] blokA_m
 *          Index of the first off-diagonal block in cblk, that is used for A.
 *
 * @param[in] blokB_m
 *          Index of the first off-diagonal block in cblk, that is used for B.
 *
 * @param[in] lrA
 *          Pointer to the low-rank representation of the column block A.
 *          Must be followed by the low-rank representation of the following blocks.
 *
 * @param[in] lrB
 *          Pointer to the low-rank representation of the column block B.
 *          Must be followed by the low-rank representation of the following blocks.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpublok_zadd_lrlr( pastix_int_t            alpha,
                   const SolverCblk       *cblkA,
                   SolverCblk             *cblkB,
                   pastix_int_t            blokA_m,
                   pastix_int_t            blokB_m,
                   const pastix_lrblock_t *lrA,
                   pastix_lrblock_t       *lrB,
                   pastix_complex64_t     *work,
                   pastix_int_t            lwork,
                   const pastix_lr_t      *lowrank )
{
    const SolverBlok   *blokA  = cblkA->fblokptr + blokA_m;
    const SolverBlok   *blokB  = cblkB->fblokptr + blokB_m;
    const SolverBlok   *lblokA = cblkA[1].fblokptr;
    const SolverBlok   *lblokB = cblkB[1].fblokptr;
    pastix_fixdbl_t     flops = 0.;
    core_zlrmm_t        params;

    assert( (cblkA->cblktype & CBLK_COMPRESSED) );
    assert( (cblkB->cblktype & CBLK_COMPRESSED) );

    params.lowrank = lowrank;
    params.transA  = PastixNoTrans; /* Unused */
    params.transB  = PastixNoTrans; /* Unused */
    params.K       = -1;            /* Unused */
    params.alpha   = alpha;
    params.A       = NULL;          /* Unused */
    params.B       = NULL;          /* Unused */
    params.beta    = 1.0;
    params.work    = work;
    params.lwork   = lwork;
    params.lwused  = 0;
    params.lock    = &(cblkB->lock);

    /* Dimensions on N */
    params.N    = cblk_colnbr( cblkA );
    params.Cn   = cblk_colnbr( cblkB );
    params.offy = cblkA->fcolnum - cblkB->fcolnum;

    do {
        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++; lrB++;
        }
        assert( blokA->fcblknm == blokB->fcblknm );
        assert( is_block_inside_fblock( blokA, blokB ) );
        assert( blokB <= lblokB );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = lrB;
        flops += core_zlradd( &params, lrA, PastixNoTrans, PASTIX_LRM3_ORTHOU );
        blokA++;
        lrA++;
    }
    while ( ( blokA < lblokA ) &&
            ( blokA[-1].fcblknm == blokA[0].fcblknm ) &&
            ( blokA[-1].lcblknm == blokA[0].lcblknm ) );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 * @param[in] blokA_m
 *          Index of the first off-diagonal block in cblk, that is used for A.
 *
 * @param[in] blokB_m
 *          Index of the first off-diagonal block in cblk, that is used for B.
 *
 * @param[inout] A
 *          The pointer to the coeftab of the blok.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          blok.ucoeftab otherwise. Must be of size blok.stride -by- blok.width
 *
 * @param[inout] B
 *          The pointer to the coeftab of the blok.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; blok.ucoeftab otherwise. Must be of size
 *          blok.stride -by- blok.width
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpublok_zadd_frfr( pastix_int_t              alpha,
                   const SolverCblk         *cblkA,
                   SolverCblk               *cblkB,
                   pastix_int_t              blokA_m,
                   pastix_int_t              blokB_m,
                   const pastix_complex64_t *A,
                   pastix_complex64_t       *B )
{
    const SolverBlok         *blokA  = cblkA->fblokptr + blokA_m;
    const SolverBlok         *blokB  = cblkB->fblokptr + blokB_m;
    const SolverBlok         *lblokA = cblkA[1].fblokptr;
    const SolverBlok         *lblokB = cblkB[1].fblokptr;
    const pastix_complex64_t *bA     = A;
    pastix_complex64_t       *bB     = B;
    pastix_int_t              lda;
    pastix_int_t              ldb;
    pastix_int_t              m;
    pastix_int_t              n     = cblk_colnbr( cblkA );
    pastix_fixdbl_t           flops = 0.;
    size_t offsetA, offsetB;

    assert( !(cblkA->cblktype & CBLK_COMPRESSED) );
    assert( !(cblkB->cblktype & CBLK_COMPRESSED) );

    /* Both cblk A and B must be stored in 2D */
    assert( cblkA->cblktype & CBLK_LAYOUT_2D );
    assert( cblkB->cblktype & CBLK_LAYOUT_2D );

    offsetA = blokA->coefind;
    offsetB = blokB->coefind;

    /* Adds blocks facing the same cblk */
    do {
        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++;
        }
        assert( blokA->fcblknm == blokB->fcblknm );
        assert( is_block_inside_fblock( blokA, blokB ) );
        assert( blokB <= lblokB );

        bA = A + blokA->coefind - offsetA;
        bB = B + blokB->coefind - offsetB;
        lda = blok_rownbr( blokA );
        ldb = blok_rownbr( blokB );

        bB = bB + ldb * ( cblkA->fcolnum - cblkB->fcolnum ) + ( blokA->frownum - blokB->frownum );
        m = lda;

        pastix_cblk_lock( cblkB );
        core_zgeadd( PastixNoTrans, m, n,
                     alpha, bA, lda,
                        1., bB, ldb );
        pastix_cblk_unlock( cblkB );
        flops += 2. * m * n;
        blokA++;
    }
    while ( ( blokA < lblokA ) &&
            ( blokA[-1].fcblknm == blokA[0].fcblknm ) &&
            ( blokA[-1].lcblknm == blokA[0].lcblknm ) );

    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in full rank format.
 *
 * The second cblk is overwritten by the sum of the two column blocks.
 *              B <- alpha * A + B
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha
 *
 * @param[in] cblkA
 *          The column block of the A matrix.
 *
 * @param[inout] cblkB
 *          The column block of the B matrix
 *          On exit, cblkB coefficient arrays are overwritten by the result of
 *          alpha * A + B.
 *
 * @param[in] blokA_m
 *          Index of the first off-diagonal block in cblk, that is used for A.
 *
 * @param[in] blokB_m
 *          Index of the first off-diagonal block in cblk, that is used for B.
 *
 * @param[inout] A
 *          The pointer to the coeftab of the blok.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          blok.ucoeftab otherwise. Must be of size blok.stride -by- blok.width
 *
 * @param[inout] B
 *          The pointer to the coeftab of the blok.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; blok.ucoeftab otherwise. Must be of size
 *          blok.stride -by- blok.width
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
pastix_fixdbl_t
cpublok_zadd( double              alpha,
              const SolverCblk   *cblkA,
              SolverCblk         *cblkB,
              pastix_int_t        blokA_m,
              pastix_int_t        blokB_m,
              const void         *A,
              void               *B,
              pastix_complex64_t *work,
              pastix_int_t        lwork,
              const pastix_lr_t  *lowrank )
{
    pastix_ktype_t ktype = PastixKernelGEADDCblkFRFR;
    pastix_fixdbl_t time, flops = 0.0;
    pastix_int_t m = cblkA->stride;
    pastix_int_t n = cblk_colnbr( cblkA );

    if ( cblkB->cblktype & CBLK_COMPRESSED ) {
        if ( cblkA->cblktype & CBLK_COMPRESSED ) {
            ktype = PastixKernelGEADDCblkLRLR;
            time  = kernel_trace_start( ktype );
            flops = cpublok_zadd_lrlr( alpha, cblkA, cblkB, blokA_m, blokB_m,
                                       A, B, work, lwork, lowrank );
        }
        else {
            ktype = PastixKernelGEADDCblkFRLR;
            time  = kernel_trace_start( ktype );
            flops = cpublok_zadd_frlr( alpha, cblkA, cblkB, blokA_m, blokB_m,
                                       A, B, work, lwork, lowrank );
        }
    }
    else {
        if ( cblkA->cblktype & CBLK_COMPRESSED ) {
            assert(0); /* We do not add a compressed cblk to a non compressed cblk */
            return 0.; /* Avoids compilation and coverity warning */
        }
        else {
            ktype = PastixKernelGEADDCblkFRFR;
            time  = kernel_trace_start( ktype );
            flops = cpublok_zadd_frfr( alpha, cblkA, cblkB, blokA_m, blokB_m, A, B );
        }
    }

    kernel_trace_stop( cblkB->fblokptr->inlast, ktype, m, n, 0, flops, time );
    return flops;
}
