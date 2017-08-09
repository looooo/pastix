/**
 *
 * @file core_zsytrfsp.c
 *
 * PaStiX kernel routines for LDL^t factorization.
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
#include "pastix_zcores.h"

#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define MAXSIZEOFBLOCKS 64
static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup kernel_blas_lapack_null
 *
 * @brief Compute the sequential static pivoting factorization of the symmetric
 * matrix n-by-n A such that A = L * D * L^t.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LDL^t factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************/
static void
core_zsytf2sp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivot,
               double              criteria )
{
    pastix_int_t k, m;
    pastix_complex64_t *Akk = A;     /* A [k  ][k  ] */
    pastix_complex64_t *Amk = A+1;   /* A [k+1][k  ] */
    pastix_complex64_t *Akm = A+lda; /* A [k  ][k+1] */
    pastix_complex64_t  alpha;

    m = n-1;
    for (k=0; k<n; k++, m--){
        if ( cabs(*Akk) < criteria ) {
            (*Akk) = (pastix_complex64_t)criteria;
            (*nbpivot)++;
        }

        alpha = 1. / (*Akk);

        /* Transpose the column before scaling */
        cblas_zcopy( m, Amk, 1, Akm, lda );

        /* Scale the diagonal to compute L((k+1):n,k) */
        cblas_zscal(m, CBLAS_SADDR( alpha ), Amk, 1 );

        alpha = -(*Akk);

        /* Move to next Akk */
        Akk += (lda+1);

        cblas_zsyrk(CblasColMajor, CblasLower, CblasNoTrans,
                    m, 1,
                    CBLAS_SADDR( alpha ), Amk, lda,
                    CBLAS_SADDR( zone  ), Akk, lda);

        /* Move to next Amk */
        Amk = Akk+1;
        Akm = Akk+lda;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the block static pivoting factorization of the symmetric
 * matrix n-by-n A such that A = L * D * L^t.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LDL^t factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[inout] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************/
void
core_zsytrfsp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivot,
               double              criteria )
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, col;
    pastix_complex64_t *Akk, *Amk, *Akm, *Amm;
    pastix_complex64_t alpha;

    /* diagonal supernode is divided into MAXSIZEOFBLOCK-by-MAXSIZEOFBLOCKS blocks */
    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    for (k=0; k<blocknbr; k++) {

        blocksize = pastix_imin(MAXSIZEOFBLOCKS, n-k*MAXSIZEOFBLOCKS);
        Akk = A+(k*MAXSIZEOFBLOCKS)*(lda+1); /* Lk,  k   */
        Amk = Akk + blocksize;               /* Lk+1,k   */
        Akm = Akk + blocksize * lda;         /* Lk,  k+1 */
        Amm = Amk + blocksize * lda;         /* Lk+1,k+1 */

        /* Factorize the diagonal block Akk*/
        core_zsytf2sp(blocksize, Akk, lda, nbpivot, criteria);

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
                        CblasTrans, CblasUnit,
                        matrixsize, blocksize,
                        CBLAS_SADDR(zone), Akk, lda,
                                           Amk, lda);

            /* Compute L(k+1:n,k) = A(k+1:n,k)D(k,k)^{-1}     */
            for(col = 0; col < blocksize; col++) {
                /* copy L(k+1+col:n,k+col)*D(k+col,k+col) into work(:,col) */
                cblas_zcopy(matrixsize, Amk + col*lda, 1,
                                        Akm + col,     lda);

                /* compute L(k+1+col:n,k+col) = A(k+1+col:n,k+col)D(k+col,k+col)^{-1} */
                alpha = 1. / *(Akk + col*(lda+1));
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
}

/**
 *******************************************************************************
 *
 * @brief Computes the LDL^t factorization of the diagonal block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting performed during the diagonal block
 *         factorization.
 *
 *******************************************************************************/
int
cpucblk_zsytrfsp1d_sytrf( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          double              criteria )
{
    pastix_int_t  ncols, stride;
    pastix_int_t  nbpivot = 0;

    ncols  = cblk->lcolnum - cblk->fcolnum + 1;
    stride = (cblk->cblktype & CBLK_LAYOUT_2D) ? ncols : cblk->stride;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        assert( cblk->fblokptr->LRblock[0].rk == -1 );
        L = cblk->fblokptr->LRblock[0].u;
        stride = ncols;

        assert( stride == cblk->fblokptr->LRblock[0].rkmax );
    }

    /*
     * Factorize diagonal block in L D L^t
     *
     *  - lower part holds L
     *  - diagonal holds D
     *  - uppert part holds (DL^t)
     */
    core_zsytrfsp( ncols, L, stride, &nbpivot, criteria );

    return nbpivot;
}

/**
 *******************************************************************************
 *
 * core_zsytrfsp1d_trsm - Apply all the trsm updates to one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit.
 *
 *******************************************************************************/
int core_zsytrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L )
{
    const SolverBlok *blok, *lblk;

    blok = cblk->fblokptr + 1; /* Firt off-diagonal block */
    lblk = cblk[1].fblokptr;   /* Next diagonal block     */

    /* if there are off-diagonal supernodes in the column */
    if ( blok < lblk )
    {
        const pastix_complex64_t *A;
        pastix_complex64_t *B;
        pastix_int_t M, N, lda, ldb, j;

        N = cblk_colnbr( cblk );
        A = L;

        if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
            lda = N;

            for(; blok < lblk; blok++) {
                M   = blok_rownbr( blok );
                B   = L + blok->coefind;
                ldb = M;

                /* Three terms version, no need to keep L and L*D */
                cblas_ztrsm( CblasColMajor, CblasRight, CblasLower,
                             CblasTrans, CblasUnit, M, N,
                             CBLAS_SADDR(zone), A, lda,
                                                B, ldb);

                for (j=0; j<N; j++)
                {
                    pastix_complex64_t alpha;
                    alpha = 1. / A[j + j * lda];
                    cblas_zscal(M, CBLAS_SADDR(alpha), B + j * ldb, 1);
                }
            }

        }
        else {
            lda = cblk->stride;
            M   = cblk->stride - N;
            B   = L + blok->coefind;
            ldb = cblk->stride;

            /* Three terms version, no need to keep L and L*D */
            cblas_ztrsm( CblasColMajor, CblasRight, CblasLower,
                         CblasTrans, CblasUnit, M, N,
                         CBLAS_SADDR(zone), A, lda,
                                            B, ldb);

            for (j=0; j<N; j++)
            {
                pastix_complex64_t alpha;
                alpha = 1. / A[j + j * lda];
                cblas_zscal(M, CBLAS_SADDR(alpha), B + j * ldb, 1);
            }
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * core_zsytrfsp1d_gemm - Computes the LDL^t factorization of one panel and
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
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
void core_zsytrfsp1d_gemm( const SolverCblk         *cblk,
                           const SolverBlok         *blok,
                                 SolverCblk         *fcblk,
                           const pastix_complex64_t *L,
                                 pastix_complex64_t *C,
                                 pastix_complex64_t *work )
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
                               -1., blokA, lda,
                                    blokB, ldb,
                                1., blokC, ldc,
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
 * @brief Compute the LDL^t factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] DLt
 *          The pointer to the upper matrix storing the coefficients the
 *          temporary DL^t product. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 * @param[in] lowrank
 *          The structure with low-rank parameters.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zsytrfsp1d_panel( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *DLt,
                          double              criteria,
                          const pastix_lr_t  *lowrank )
{
    pastix_int_t nbpivot;
    (void)lowrank;

    nbpivot = cpucblk_zsytrfsp1d_sytrf( cblk, L, criteria );

    /*
     * We exploit the fact that (DL^t) is stored in the upper triangle part of L
     */
    cpucblk_ztrsmsp( PastixLCoef, PastixRight, PastixUpper,
                     PastixNoTrans, PastixNonUnit,
                     cblk, L, L, lowrank );

    if ( DLt != NULL ) {
        /* Copy L into the temporary buffer and multiply by D */
        cpucblk_zscalo( PastixNoTrans, cblk, DLt );
    }
    return nbpivot;
}


/**
 *******************************************************************************
 *
 * @brief Perform the LDL^t factorization of a given panel and apply all its
 * updates.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 * @param[in] DLt
 *          Temporary memory buffer to store the transpose of DLt.
 *
 * @param[in] work2
 *          Temporary memory buffer for U factors.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zsytrfsp1d( SolverMatrix       *solvmtx,
                    SolverCblk         *cblk,
                    double              criteria,
                    pastix_complex64_t *DLt,
                    pastix_complex64_t *work2 )
{
    pastix_complex64_t *L = cblk->lcoeftab;
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivot;

    /* if there are off-diagonal supernodes in the column */
    nbpivot = cpucblk_zsytrfsp1d_panel(cblk, L, DLt, criteria, &solvmtx->lowrank );

    blok = cblk->fblokptr+1;   /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        /* Update on L */
        if (DLt == NULL) {
            core_zsytrfsp1d_gemm( cblk, blok, fcblk,
                                  L, fcblk->lcoeftab,
                                  work2 );
        }
        else {
            cpucblk_zgemmsp( PastixLCoef, PastixUCoef, PastixTrans,
                             cblk, blok, fcblk,
                             L, DLt, fcblk->lcoeftab,
                             work2, &solvmtx->lowrank );
        }
        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }

    return nbpivot;
}
