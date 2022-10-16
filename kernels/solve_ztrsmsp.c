/**
 *
 * @file solve_ztrsmsp.c
 *
 * PaStiX solve kernels routines
 *
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2021-03-05
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
static pastix_complex64_t zzero =  0.0;
static pastix_complex64_t  zone =  1.0;
static pastix_complex64_t mzone = -1.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Apply a solve trsm update related to a diagonal block of the matrix A.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify the side parameter of the TRSM.
 *
 * @param[in] uplo
 *          Specify the uplo parameter of the TRSM.
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] nrhs
 *          The number of right hand side.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the A and B matrix.
 *
 * @param[in] dataA
 *          The pointer to the correct representation of the data of A.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[inout] B
 *          The pointer to the matrix B, that is a portion of the right hand
 *          side to solve.
 *
 * @param[in] ldb
 *          The leading dimension of B.
 *
 *******************************************************************************/
void
solve_blok_ztrsm( pastix_side_t       side,
                  pastix_uplo_t       uplo,
                  pastix_trans_t      trans,
                  pastix_diag_t       diag,
                  const SolverCblk   *cblk,
                  int                 nrhs,
                  const void         *dataA,
                  pastix_complex64_t *b,
                  int                 ldb )
{
    pastix_int_t        n, lda;
    pastix_lrblock_t   *lrA;
    pastix_complex64_t *A;

    n = cblk_colnbr( cblk );

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        lrA = (pastix_lrblock_t *)dataA;
        assert( lrA->rk == -1 );
        A   = lrA->u;
        lda = n;
    }
    else {
        A   = (pastix_complex64_t *)dataA;
        lda = (cblk->cblktype & CBLK_LAYOUT_2D) ? n : cblk->stride;
    }

    cblas_ztrsm(
        CblasColMajor, (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)trans, (CBLAS_DIAG)diag,
        n, nrhs,
        CBLAS_SADDR(zone), A, lda,
                           b, ldb );
}

/**
 *******************************************************************************
 *
 * @brief Apply a solve gemm update related to a single block of the matrix A.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the blok parameter belongs to cblk (PastixLeft), or
 *          to fcbk (PastixRight).
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] nrhs
 *          The number of right hand side.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the B matrix.
 *
 * @param[in] blok
 *          The blok structure that corresponds to the A matrix, and that
 *          belongs either to cblk or fcbk depending on the side parameter.
 *
 * @param[inout] fcbk
 *          The cblk structure that corresponds to the C matrix.
 *
 * @param[in] dataA
 *          The pointer to the correct representation of the data of A.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[in] B
 *          The pointer to the matrix B, that is a portion of the right hand
 *          side.
 *
 * @param[in] ldb
 *          The leading dimension of B.
 *
 * @param[inout] C
 *          The pointer to the matrix C, that is the updated portion of the
 *          right hand side.
 *
 * @param[in] ldb
 *          The leading dimension of C.
 *
 *******************************************************************************/
void
solve_blok_zgemm( pastix_side_t             side,
                  pastix_trans_t            trans,
                  pastix_int_t              nrhs,
                  const SolverCblk         *cblk,
                  const SolverBlok         *blok,
                  SolverCblk               *fcbk,
                  const void               *dataA,
                  const pastix_complex64_t *B,
                  pastix_int_t              ldb,
                  pastix_complex64_t       *C,
                  pastix_int_t              ldc )
{
    pastix_int_t      m, n, lda;
    pastix_int_t      offB, offC;
    const SolverCblk *bowner;

    if ( side == PastixLeft ) {
        /*
         * Blok should belong to cblk
         */
        bowner = cblk;

        m = blok_rownbr( blok );
        n = cblk_colnbr( cblk );
        lda = m;

        offB = 0;
        offC = blok->frownum - fcbk->fcolnum;
    }
    else {
        /*
         * Blok should belong to fcbk
         */
        bowner = fcbk;

        m = cblk_colnbr( fcbk );
        n = blok_rownbr( blok );
        lda = n;

        offB = blok->frownum - cblk->fcolnum;
        offC = 0;
    }

    assert( (blok > bowner[0].fblokptr) &&
            (blok < bowner[1].fblokptr) );

    if ( bowner->cblktype & CBLK_COMPRESSED ) {
        const pastix_lrblock_t *lrA = dataA;
        pastix_complex64_t     *tmp;

        switch (lrA->rk){
        case 0:
            break;
        case -1:
            pastix_cblk_lock( fcbk );
            cblas_zgemm(
                CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                m, nrhs, n,
                CBLAS_SADDR(mzone), lrA->u,   lda,
                                    B + offB, ldb,
                CBLAS_SADDR(zone),  C + offC, ldc );
            pastix_cblk_unlock( fcbk );
            break;
        default:
            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);
            if (trans == PastixNoTrans) {
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    lrA->rk, nrhs, n,
                    CBLAS_SADDR(zone),  lrA->v,   lrA->rkmax,
                                        B + offB, ldb,
                    CBLAS_SADDR(zzero), tmp,      lrA->rk );

                pastix_cblk_lock( fcbk );
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    m, nrhs, lrA->rk,
                    CBLAS_SADDR(mzone), lrA->u,   lda,
                                        tmp,      lrA->rk,
                    CBLAS_SADDR(zone),  C + offC, ldc );
                pastix_cblk_unlock( fcbk );
            }
            else {
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    lrA->rk, nrhs, n,
                    CBLAS_SADDR(zone),  lrA->u,   lda,
                                        B + offB, ldb,
                    CBLAS_SADDR(zzero), tmp,      lrA->rk );

                pastix_cblk_lock( fcbk );
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    m, nrhs, lrA->rk,
                    CBLAS_SADDR(mzone), lrA->v,   lrA->rkmax,
                                        tmp,      lrA->rk,
                    CBLAS_SADDR(zone),  C + offC, ldc );
                pastix_cblk_unlock( fcbk );
            }
            memFree_null(tmp);
            break;
        }
    }
    else{
        const pastix_complex64_t *A = dataA;
        lda = (bowner->cblktype & CBLK_LAYOUT_2D) ? lda : bowner->stride;

        pastix_cblk_lock( fcbk );
        cblas_zgemm(
            CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
            m, nrhs, n,
            CBLAS_SADDR(mzone), A,        lda,
                                B + offB, ldb,
            CBLAS_SADDR(zone),  C + offC, ldc );
        pastix_cblk_unlock( fcbk );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply a forward solve related to one cblk to all the right hand side.
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
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
solve_cblk_ztrsmsp_forward( pastix_solv_mode_t  mode,
                            pastix_side_t       side,
                            pastix_uplo_t       uplo,
                            pastix_trans_t      trans,
                            pastix_diag_t       diag,
                            const SolverMatrix *datacode,
                            const SolverCblk   *cblk,
                            pastix_rhs_t        rhsb )
{
    SolverCblk       *fcbk;
    const SolverBlok *blok;
    pastix_trans_t    tA;
    pastix_coefside_t cs;
    const void               *dataA = NULL;
    const pastix_lrblock_t   *lrA;
    const pastix_complex64_t *A;
    pastix_complex64_t       *B, *C;
    pastix_int_t              ldb, ldc;

    if ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixTrans;
        cs = PastixUCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else if ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else if ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixNoTrans;
        cs = PastixUCoef;

        /* We do not handle conjtrans in complex as we store U^t */
#if defined(PRECISION_z) || defined(PRECISION_c)
        assert( trans != PastixConjTrans );
#endif
    }
    else if ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;
    }
    else {
        /* This correspond to case treated in backward trsm */
        assert(0);
        return;
    }

    assert( !( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) );

    if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) ) {
        return;
    }

    B   = rhsb->b;
    B   = B + cblk->lcolidx;
    ldb = rhsb->ld;

    /* Solve the diagonal block */
    solve_blok_ztrsm( side, PastixLower,
                      tA, diag, cblk, rhsb->n,
                      cblk_getdata( cblk, cs ),
                      B, ldb );

    /* Apply the update */
    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
        fcbk = datacode->cblktab + blok->fcblknm;

        if ( (fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeLocal) ) {
            return;
        }
        assert( !(fcbk->cblktype & CBLK_RECV) );

        /*
         * Make sure we get the correct pointer to the lrA, or to the right position in [lu]coeftab
         */
        dataA = cblk_getdata( cblk, cs );
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            lrA = dataA;
            lrA += (blok - cblk->fblokptr);
            dataA = lrA;
        }
        else {
            A = dataA;
            A += blok->coefind;
            dataA = A;
        }

        /*
         * Make sure we get the correct pointer for the C matrix.
         */
        if ( fcbk->cblktype & CBLK_FANIN ) {
            C   = rhsb->cblkb[ - fcbk->bcscnum - 1 ];
            ldc = cblk_colnbr( fcbk );
            if ( C == NULL ) {
                rhsb->cblkb[ - fcbk->bcscnum - 1 ] = calloc( ldc * rhsb->n, sizeof( pastix_complex64_t ) );
                C = rhsb->cblkb[ - fcbk->bcscnum - 1 ];
            }
        }
        else {
            C   = rhsb->b;
            C   = C + fcbk->lcolidx;
            ldc = rhsb->ld;
        }

        solve_blok_zgemm( PastixLeft, tA, rhsb->n,
                          cblk, blok, fcbk,
                          dataA, B, ldb, C, ldc );
        pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply a backward solve related to one cblk to all the right hand side.
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
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
solve_cblk_ztrsmsp_backward( pastix_solv_mode_t  mode,
                             pastix_side_t       side,
                             pastix_uplo_t       uplo,
                             pastix_trans_t      trans,
                             pastix_diag_t       diag,
                             const SolverMatrix *datacode,
                             SolverCblk         *cblk,
                             pastix_rhs_t        rhsb )
{
    SolverCblk       *fcbk;
    const SolverBlok *blok;
    pastix_int_t      j;
    pastix_trans_t    tA;
    pastix_coefside_t cs;
    const void               *dataA = NULL;
    const pastix_lrblock_t   *lrA;
    const pastix_complex64_t *A;
    pastix_complex64_t       *B, *C;
    pastix_int_t              ldb, ldc;

    /*
     *  Left / Upper / NoTrans (Backward)
     */
    if ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixTrans;
        cs = PastixUCoef;
    }
    else if ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;
    }
    else if ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) {
        /*  We store U^t, so we swap uplo and trans */
        tA = PastixNoTrans;
        cs = PastixUCoef;

        /* Right is not handled yet */
        assert( 0 );

        /* We do not handle conjtrans in complex as we store U^t */
        assert( trans != PastixConjTrans );
    }
    else if ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) {
        tA = trans;
        cs = PastixLCoef;

        /* Right is not handled yet */
        assert( 0 );
    }
    else {
        /* This correspond to case treated in forward trsm */
        assert(0);
        return;
    }

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

    /*
     * Make sure we get the correct pointer for the B matrix.
     */
    assert( !(cblk->cblktype & CBLK_RECV) );
    if ( cblk->cblktype & CBLK_FANIN ) {
        B   = rhsb->cblkb[ - cblk->bcscnum - 1 ];
        ldb = cblk_colnbr( cblk );
    }
    else {
        B   = rhsb->b;
        B   = B + cblk->lcolidx;
        ldb = rhsb->ld;
    }

    if ( !(cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) &&
         (!(cblk->cblktype & CBLK_IN_SCHUR) || (mode == PastixSolvModeSchur)) )
    {
        /* Solve the diagonal block */
        solve_blok_ztrsm( side, PastixLower, tA, diag,
                          cblk, rhsb->n,
                          cblk_getdata( cblk, cs ),
                          B, ldb );

    }

    /* Apply the update */
    for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
        blok = datacode->bloktab + datacode->browtab[j];
        fcbk = datacode->cblktab + blok->lcblknm;

        if ( ((fcbk->cblktype & CBLK_IN_SCHUR) && (mode == PastixSolvModeInterface))
           || (fcbk->cblktype & CBLK_RECV) ) {
            continue;
        }
        assert( !(fcbk->cblktype & CBLK_FANIN) );

        /*
         * Make sure we get the correct pointer to the lrA, or to the right position in [lu]coeftab
         */
        dataA = cblk_getdata( fcbk, cs );
        if ( fcbk->cblktype & CBLK_COMPRESSED ) {
            lrA = dataA;
            lrA += (blok - fcbk->fblokptr);
            dataA = lrA;
        }
        else {
            A = dataA;
            A += blok->coefind;
            dataA = A;
        }

        /*
         * Make sure we get the correct pointer for the C matrix.
         */
        C   = rhsb->b;
        C   = C + fcbk->lcolidx;
        ldc = rhsb->ld;

        solve_blok_zgemm( PastixRight, tA, rhsb->n,
                          cblk, blok, fcbk,
                          dataA, B, ldb, C, ldc );
        pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
    }

    if ( cblk->cblktype & CBLK_FANIN ) {
        memFree_null( rhsb->cblkb[ - cblk->bcscnum - 1 ] );
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
solve_cblk_zdiag( const SolverCblk   *cblk,
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
        A = (pastix_complex64_t*)(cblk->fblokptr->LRblock[0]->u);
        assert( cblk->fblokptr->LRblock[0]->rkmax == lda );
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
