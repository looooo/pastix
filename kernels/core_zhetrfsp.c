/**
 *
 * @file core_zhetrfsp.c
 *
 * PaStiX kernel routines for LDL^h factorization.
 *
 * @copyright 2011-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Alycia Lisito
 * @author Nolan Bredel
 * @date 2023-11-07
 * @precisions normal z -> c
 *
 **/
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "kernels_trace.h"

#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define MAXSIZEOFBLOCKS 64
static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t mzone = -1.0;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup kernel_blas_lapack_null
 *
 * @brief Compute the sequential static pivoting factorization of the hermitian
 * matrix n-by-n A such that A = L * D * conj(L^t).
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LDL^h factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivots
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criterion
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************/
static inline void
core_zhetf2sp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivots,
               double              criterion )
{
    pastix_int_t k, m, ret;
    pastix_complex64_t *Akk = A;     /* A [k  ][k  ] */
    pastix_complex64_t *Amk = A+1;   /* A [k+1][k  ] */
    pastix_complex64_t *Akm = A+lda; /* A [k  ][k+1] */
    pastix_complex64_t  zalpha;
    double dalpha;

    m = n-1;
    for (k=0; k<n; k++, m--){
        if ( cabs(*Akk) < criterion ) {
            if ( creal(*Akk) < 0. ) {
                *Akk = (pastix_complex64_t)(-criterion);
            }
            else {
                *Akk = (pastix_complex64_t)criterion;
            }
            (*nbpivots)++;
        }

        zalpha = 1.0 / (*Akk);

        /* Transpose the column before scaling */
        cblas_zcopy( m, Amk, 1, Akm, lda );
        ret = LAPACKE_zlacgv_work( m, Akm, 1 );
        assert( ret == 0 );

        /* Scale the diagonal to compute L((k+1):n,k) */
        cblas_zscal(m, CBLAS_SADDR( zalpha ), Amk, 1 );

        dalpha = -1.0 * creal(*Akk);

        /* Move to next Akk */
        Akk += (lda+1);

        cblas_zher(CblasColMajor, CblasLower,
                   m, dalpha,
                   Amk, 1,
                   Akk, lda);

        /* Move to next Amk */
        Amk = Akk+1;
        Akm = Akk+lda;
    }
    (void)ret;
}

/**
 *******************************************************************************
 *
 * @brief Compute the block static pivoting factorization of the hermitian
 * matrix n-by-n A such that A = L * D * conj(L^t).
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LDL^h factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivots
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criterion
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************/
void
core_zhetrfsp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivots,
               double              criterion )
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, col, ret;
    pastix_complex64_t *Akk, *Amk, *Akm, *Amm;
    pastix_complex64_t alpha;

    /* diagonal supernode is divided into MAXSIZEOFBLOCK-by-MAXSIZEOFBLOCKS blocks */
    blocknbr = pastix_iceil( n, MAXSIZEOFBLOCKS );

    for (k=0; k<blocknbr; k++) {

        blocksize = pastix_imin(MAXSIZEOFBLOCKS, n-k*MAXSIZEOFBLOCKS);
        Akk = A+(k*MAXSIZEOFBLOCKS)*(lda+1); /* Lk,  k   */
        Amk = Akk + blocksize;               /* Lk+1,k   */
        Akm = Akk + blocksize * lda;         /* Lk,  k+1 */
        Amm = Amk + blocksize * lda;         /* Lk+1,k+1 */

        /* Factorize the diagonal block Akk*/
        core_zhetf2sp(blocksize, Akk, lda, nbpivots, criterion);

        if ((k*MAXSIZEOFBLOCKS+blocksize) < n) {

            matrixsize = n-(k*MAXSIZEOFBLOCKS+blocksize);

            /*
             * Solve the lower rectangle below the diagonal block
             *      L(k+1:n,k) = (L(k,k) D(k,k))^{-1} A(k+1:n,k)
             */
            /* 1) Compute A(k+1:n,k) = A(k+1:n,k)L(k,k)^{-T} = D(k,k)L(k+1:n,k) */
                        /* input: L(k,k) in tmp, A(k+1:n,k) in tmp1   */
                        /* output: A(k+1:n,k) in tmp1                 */
            cblas_ztrsm(CblasColMajor,
                        CblasRight, CblasLower,
                        CblasConjTrans, CblasUnit,
                        matrixsize, blocksize,
                        CBLAS_SADDR(zone), Akk, lda,
                                           Amk, lda);

            /* Compute L(k+1:n,k) = A(k+1:n,k)D(k,k)^{-1}     */
            for(col = 0; col < blocksize; col++) {
                /* copy L(k+1+col:n,k+col)*D(k+col,k+col) into work(:,col) */
                cblas_zcopy(matrixsize, Amk + col*lda, 1,
                                        Akm + col,     lda);
                ret = LAPACKE_zlacgv_work( matrixsize, Akm + col, lda );
                assert( ret == 0 );

                /* compute L(k+1+col:n,k+col) = A(k+1+col:n,k+col)D(k+col,k+col)^{-1} */
                alpha = 1.0 / *(Akk + col*(lda+1));
                cblas_zscal( matrixsize, CBLAS_SADDR(alpha),
                             Amk + col*lda, 1 );
            }

            /* Update A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - (L(k+1:n,k)*D(k,k))*L(k+1:n,k)^T */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), Amk, lda,
                                            Akm, lda,
                        CBLAS_SADDR(zone),  Amm, lda);
        }
    }
    (void)ret;
}

/**
 *******************************************************************************
 *
 * @brief Computes the LDL^h factorization of the diagonal block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] dataL
 *          The pointer to the correct representation of lower part of the data.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting performed during the diagonal block
 *         factorization.
 *
 *******************************************************************************/
int
cpucblk_zhetrfsp1d_hetrf( SolverMatrix *solvmtx,
                          SolverCblk   *cblk,
                          void         *dataL )
{
    pastix_int_t  ncols, stride;
    pastix_int_t  nbpivots = 0;
    pastix_fixdbl_t time, flops;
    pastix_complex64_t *L;
    pastix_lrblock_t   *lrL;
    double criterion = solvmtx->diagthreshold;

    time = kernel_trace_start( PastixKernelHETRF );

    ncols  = cblk->lcolnum - cblk->fcolnum + 1;
    stride = (cblk->cblktype & CBLK_LAYOUT_2D) ? ncols : cblk->stride;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        /* dataL is a LRblock */
        lrL = (pastix_lrblock_t *)dataL;
        L   = lrL->u;
        stride = ncols;

        assert( lrL->rk == -1 );
        assert( stride == lrL->rkmax );
    } else {
        L = (pastix_complex64_t *)dataL;
    }

    /*
     * Factorize diagonal block in L D L^h
     *
     *  - lower part holds L
     *  - diagonal holds D
     *  - uppert part holds (DL^h)
     */
    flops = FLOPS_ZHETRF( ncols );
    kernel_trace_start_lvl2( PastixKernelLvl2HETRF );
    core_zhetrfsp( ncols, L, stride, &nbpivots, criterion );
    kernel_trace_stop_lvl2( flops );

    kernel_trace_stop( cblk->fblokptr->inlast, PastixKernelHETRF, ncols, 0, 0, flops, time );

    if ( nbpivots ) {
        pastix_atomic_add_32b( &(solvmtx->nbpivots), nbpivots );
    }
    return nbpivots;
}

/**
 *******************************************************************************
 *
 * core_zhetrfsp1d_gemm - Computes the LDL^h factorization of one panel and
 * apply all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] blok
 *          The pointer to the data structure that describes the blok from which
 *          we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on
 *          which we compute the contributions. Next column blok must be
 *          accessible through fcblk[1].
 *
 * @param[inout] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] C
 *          The pointer to the matrix storing the coefficients of the
 *          target.
 *
 * @param[inout] work
 *          Temporary buffer used in core_zgemdm().
 *
 *******************************************************************************/
void core_zhetrfsp1d_gemm( const SolverCblk         *cblk,
                           const SolverBlok         *blok,
                           SolverCblk               *fcblk,
                           const pastix_complex64_t *L,
                           pastix_complex64_t       *C,
                           pastix_complex64_t       *work )
{
    const SolverBlok *iterblok;
    const SolverBlok *fblok;
    const SolverBlok *lblok;
    const pastix_complex64_t *blokA;
    const pastix_complex64_t *blokB;
    const pastix_complex64_t *blokD;
    pastix_complex64_t *blokC;

    pastix_int_t M, N, K, lda, ldb, ldc, ldd;

    /* Get the panel update dimensions */
    K = cblk_colnbr( cblk );
    N = blok_rownbr( blok );

    /* Get info for diagonal, and the B block */
    blokD = L;
    blokB = L + blok->coefind;
    if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
        ldb = N;
        ldd = K + 1;
    }
    else {
        ldb = cblk->stride;
        ldd = cblk->stride + 1;
    }

    /*
     * Add contribution to C in fcblk:
     *    Get the first facing block of the distant panel, and the last block of
     *    the current cblk
     */
    fblok = fcblk->fblokptr;
    lblok = cblk[1].fblokptr;

    for (iterblok=blok; iterblok<lblok; iterblok++) {

        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        /* Get the A block and its dimensions */
        M     = blok_rownbr( iterblok );
        blokA = L + iterblok->coefind;
        lda   = (cblk->cblktype & CBLK_LAYOUT_2D) ? M : cblk->stride;

        /* Get the C block */
        ldc = (fcblk->cblktype & CBLK_LAYOUT_2D) ? blok_rownbr(fblok) : fcblk->stride;

        blokC = C + fblok->coefind
            + iterblok->frownum - fblok->frownum
            + (blok->frownum - fcblk->fcolnum) * ldc;

        {
            pastix_int_t ldw;
            int ret;

            /* Compute ldw which should never be larger than SOLVE_COEFMAX */
            ldw = (M+1) * K;

            pastix_cblk_lock( fcblk );
            ret = core_zgemdm( PastixNoTrans, PastixConjTrans,
                               M, N, K,
                               -1.0, blokA, lda,
                                     blokB, ldb,
                                1.0, blokC, ldc,
                                     blokD, ldd,
                               work, ldw );
            pastix_cblk_unlock( fcblk );
            assert(ret == PASTIX_SUCCESS);
            (void)ret;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the LDL^h factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *           The pointer to the correct representation of lower part of the data.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[inout] DLh
 *           The pointer to the correct representation of DLh matrix
 *           (stored in the upper part by default).
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zhetrfsp1d_panel( SolverMatrix *solvmtx,
                          SolverCblk   *cblk,
                          void         *L,
                          void         *DLh )
{
    pastix_int_t nbpivots;
    nbpivots = cpucblk_zhetrfsp1d_hetrf( solvmtx, cblk, L );

    /*
     * We exploit the fact that (DL^h) is stored in the upper triangle part of L
     */
    cpucblk_ztrsmsp( PastixRight, PastixUpper,
                     PastixNoTrans, PastixNonUnit,
                     cblk, L, L, &(solvmtx->lowrank) );

    if ( (DLh != NULL) && (cblk->cblktype & CBLK_LAYOUT_2D) ) {

        /* Copy L into the temporary buffer and multiply by D */
        cpucblk_zscalo( PastixConjTrans, cblk, L, DLh );
    }
    return nbpivots;
}


/**
 *******************************************************************************
 *
 * @brief Perform the LDL^h factorization of a given panel and apply all its
 * updates.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in] DLh
 *          Temporary memory buffer to store the conjugate transpose of DLh.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zhetrfsp1d( SolverMatrix       *solvmtx,
                    SolverCblk         *cblk,
                    pastix_complex64_t *DLh,
                    pastix_complex64_t *work,
                    pastix_int_t        lwork )
{
    void        *dataL   = cblk_getdataL( cblk );
    void        *dataDLh = cblk_getdataU( cblk );
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivots;

    if ( !(cblk->cblktype & CBLK_LAYOUT_2D) ) {
        DLh = NULL;
    }
    else {
        if (cblk->cblktype & CBLK_COMPRESSED) {
            cpucblk_zalloc_lrws( cblk, dataDLh, DLh );
        }
        else {
            assert( dataDLh == NULL );
            dataDLh = DLh;
        }
    }

    /* if there are off-diagonal supernodes in the column */
    nbpivots = cpucblk_zhetrfsp1d_panel( solvmtx, cblk, dataL, dataDLh );

    blok = cblk->fblokptr+1; /* this diagonal block */
    lblk = cblk[1].fblokptr; /* the next diagonal block */

    for( ; blok < lblk; blok++ )
    {
        fcblk = solvmtx->cblktab + blok->fcblknm;

        if ( fcblk->cblktype & CBLK_FANIN ) {
            cpucblk_zalloc( PastixLCoef, fcblk );
        }

        /* Update on L */
        if ( DLh == NULL ) {
            core_zhetrfsp1d_gemm( cblk, blok, fcblk,
                                  dataL, cblk_getdataL( fcblk ),
                                  work );
        }
        else {
            cpucblk_zgemmsp( PastixLCoef, PastixTrans,
                             cblk, blok, fcblk,
                             dataL, dataDLh, cblk_getdataL( fcblk ),
                             work, lwork, &(solvmtx->lowrank) );
        }
        cpucblk_zrelease_deps( PastixLCoef, solvmtx, cblk, fcblk );
    }

    return nbpivots;
}

/**
 *******************************************************************************
 *
 * @brief Perform the LDL^h factorization of a given panel and submit tasks
 * for the subsequent updates.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zhetrfsp1dplus( SolverMatrix *solvmtx,
                        SolverCblk   *cblk )
{
    void           *dataL   = cblk_getdataL( cblk );
    SolverCblk     *fcblk;
    SolverBlok     *blok, *lblk;
    pastix_int_t    i, nbpivots;
    pastix_queue_t *queue    = solvmtx->computeQueue[ cblk->threadid ];

    assert( cblk->cblktype & CBLK_TASKS_2D );
    nbpivots = cpucblk_zhetrfsp1d_panel( solvmtx, cblk, dataL, NULL );

    blok = cblk->fblokptr + 1; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    /* if there are off-diagonal supernodes in the column */
    for( i=0; blok < lblk; i++, blok++ )
    {
        fcblk = solvmtx->cblktab + blok->fcblknm;

        assert( !(fcblk->cblktype & CBLK_RECV) );
        pqueuePush1( queue, - (blok - solvmtx->bloktab) - 1, cblk->priority + i );

        /* Skip blocks facing the same cblk */
        while ( ( blok < lblk ) &&
                ( blok[0].fcblknm == blok[1].fcblknm ) &&
                ( blok[0].lcblknm == blok[1].lcblknm ) )
        {
            blok++;
        }
    }

    return nbpivots;
}

/**
 *******************************************************************************
 *
 * @brief Apply the updates of the LDL^h factorisation of a given panel.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          Solver Matrix structure of the problem
 *
 * @param[in] blok
 *          Pointer to the blok where the update start.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] lwork
 *          Temporary workspace dimension.
 *
 *******************************************************************************/
void
cpucblk_zhetrfsp1dplus_update( SolverMatrix       *solvmtx,
                               SolverBlok         *blok,
                               pastix_complex64_t *work )
{
    SolverCblk *cblk  = solvmtx->cblktab + blok->lcblknm;
    SolverCblk *fcbk  = solvmtx->cblktab + blok->fcblknm;
    SolverBlok *lblk  = cblk[1].fblokptr;   /* the next diagonal block */
    void       *dataL = cblk_getdataL( cblk );

    if ( fcbk->cblktype & CBLK_FANIN ) {
        cpucblk_zalloc( PastixLCoef, fcbk );
    }

    do
    {
        /* Update on L (3 terms) */
        core_zhetrfsp1d_gemm( cblk, blok, fcbk,
                              dataL, fcbk->lcoeftab,
                              work );

        cpucblk_zrelease_deps( PastixLCoef, solvmtx, cblk, fcbk );
        blok++;
    }
    while ( ( blok < lblk ) &&
            ( blok[-1].fcblknm == blok[0].fcblknm ) &&
            ( blok[-1].lcblknm == blok[0].lcblknm ) );
}
