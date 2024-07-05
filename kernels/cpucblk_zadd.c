/**
 *
 * @file cpucblk_zadd.c
 *
 * Precision dependent routines to add different cblks.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Nolan Bredel
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
 * @brief Add a column blok in full rank format to a column blok in low rank
 * format.
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
 * @param[inout] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
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
cpucblk_zadd_frlr( pastix_complex64_t        alpha,
                   const SolverCblk         *cblkA,
                   SolverCblk               *cblkB,
                   const pastix_complex64_t *A,
                   pastix_lrblock_t         *lrB,
                   pastix_complex64_t       *work,
                   pastix_int_t              lwork,
                   const pastix_lr_t        *lowrank )
{
    const SolverBlok   *blokA  = cblkA->fblokptr;
    const SolverBlok   *blokB  = cblkB->fblokptr;
    const SolverBlok   *lblokA = cblkA[1].fblokptr;
    const SolverBlok   *lblokB = cblkB[1].fblokptr;
    pastix_fixdbl_t     flops = 0.;
    core_zlrmm_t        params;
    pastix_lrblock_t    lrA;

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

    for (; blokA < lblokA; blokA++) {

        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++; lrB++;
        }

        assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

        lrA.u     = (pastix_complex64_t*)A + blokA->coefind;
        lrA.rkmax = blok_rownbr( blokA );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = lrB;

        flops += core_zlradd( &params, &lrA,
                              PastixNoTrans, 0 );
    }
    return flops;
}

/**
 *******************************************************************************
 *
 * @brief Add two column bloks in low rank format.
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
cpucblk_zadd_lrlr( pastix_complex64_t      alpha,
                   const SolverCblk       *cblkA,
                   SolverCblk             *cblkB,
                   const pastix_lrblock_t *lrA,
                   pastix_lrblock_t       *lrB,
                   pastix_complex64_t     *work,
                   pastix_int_t            lwork,
                   const pastix_lr_t      *lowrank )
{
    const SolverBlok   *blokA  = cblkA->fblokptr;
    const SolverBlok   *blokB  = cblkB->fblokptr;
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

    for (; blokA < lblokA; blokA++, lrA++) {

        /* Find facing bloknum */
        while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
            blokB++; lrB++;
        }

        assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

        /* Dimensions on M */
        params.M    = blok_rownbr( blokA );
        params.Cm   = blok_rownbr( blokB );
        params.offx = blokA->frownum - blokB->frownum;
        params.C    = lrB;
        flops += core_zlradd( &params, lrA, PastixNoTrans, PASTIX_LRM3_ORTHOU );
    }
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
 * @param[inout] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] B
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return The number of flops of the operation.
 *
 *******************************************************************************/
static inline pastix_fixdbl_t
cpucblk_zadd_frfr( pastix_complex64_t        alpha,
                   const SolverCblk         *cblkA,
                   SolverCblk               *cblkB,
                   const pastix_complex64_t *A,
                   pastix_complex64_t       *B )
{
    pastix_int_t    n = cblk_colnbr( cblkA );
    pastix_int_t    m = cblkA->stride;
    pastix_fixdbl_t flops = m * n;

    assert( !(cblkA->cblktype & CBLK_COMPRESSED) );
    assert( !(cblkB->cblktype & CBLK_COMPRESSED) );

    assert( (A != NULL) && (B != NULL) );

    /* If the cblk matches */
    if ( (n == cblk_colnbr( cblkB )) &&
         (m == cblkB->stride) ) {

        pastix_cblk_lock( cblkB );
        core_zgeadd( PastixNoTrans, m, n,
                     alpha, A, m,
                        1., B, m );
        pastix_cblk_unlock( cblkB );
    }
    else {
        const pastix_complex64_t *bA;
        pastix_complex64_t       *bB;
        const SolverBlok         *blokA  = cblkA->fblokptr;
        const SolverBlok         *blokB  = cblkB->fblokptr;
        const SolverBlok         *lblokA = cblkA[1].fblokptr;
        const SolverBlok         *lblokB = cblkB[1].fblokptr;
        pastix_int_t              lda, ldb;

        /* Both cblk A and B must be stored in 2D */
        assert( cblkA->cblktype & CBLK_LAYOUT_2D );
        assert( cblkB->cblktype & CBLK_LAYOUT_2D );

        for (; blokA < lblokA; blokA++) {

            /* Find facing bloknum */
            while ( !is_block_inside_fblock( blokA, blokB ) && (blokB < lblokB) ) {
                blokB++;
            }

            assert( is_block_inside_fblock( blokA, blokB ) && (blokB <= lblokB) );

            bA = A + blokA->coefind;
            bB = B + blokB->coefind;
            lda = blok_rownbr( blokA );
            ldb = blok_rownbr( blokB );

            bB = bB + ldb * ( cblkA->fcolnum - cblkB->fcolnum ) + ( blokA->frownum - blokB->frownum );
            m = lda;

            pastix_cblk_lock( cblkB );
            core_zgeadd( PastixNoTrans, m, n,
                         alpha, bA, lda,
                            1., bB, ldb );
            pastix_cblk_unlock( cblkB );
        }
    }
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
 * @param[inout] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
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
cpucblk_zadd( pastix_complex64_t  alpha,
              const SolverCblk   *cblkA,
              SolverCblk         *cblkB,
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
            flops = cpucblk_zadd_lrlr( alpha, cblkA, cblkB,
                                       A, B, work, lwork, lowrank );
        }
        else {
            ktype = PastixKernelGEADDCblkFRLR;
            time  = kernel_trace_start( ktype );
            flops = cpucblk_zadd_frlr( alpha, cblkA, cblkB,
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
            flops = cpucblk_zadd_frfr( alpha, cblkA, cblkB, A, B );
        }
    }

    kernel_trace_stop( cblkB->fblokptr->inlast, ktype, m, n, 0, flops, time );
    return flops;
}

