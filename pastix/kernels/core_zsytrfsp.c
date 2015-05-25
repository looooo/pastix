/**
 *
 * @file core_zsytrfsp.c
 *
 *  PaStiX kernel routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix_zcores.h"
#include <cblas.h>
#include "../blend/solver.h"

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytf2sp - Computes the sequential static pivoting
 * factorization of the hermitian matrix n-by-n A such that A = L * D * L^t.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with Cholesky factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[in,out] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
static void core_zsytf2sp(pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        lda,
                          pastix_int_t       *nbpivot,
                          double              criteria )
{
    pastix_int_t k;
    pastix_complex64_t *Akk = A;   /* A [k  ][k] */
    pastix_complex64_t *Amk = A+1; /* A [k+1][k] */
    pastix_complex64_t  alpha;

    for (k=0; k<n; k++){
        if ( cabs(*Akk) < criteria ) {
            (*Akk) = (pastix_complex64_t)criteria;
            (*nbpivot)++;
        }

        alpha = 1. / (*Akk);

        /* Scale the diagonal to compute L((k+1):n,k) */
        cblas_zscal(n-k-1, CBLAS_SADDR( alpha ), Amk, 1 );

        alpha = -(*Akk);

        /* Move to next Akk */
        Akk += (lda+1);

        /* TODO: check */
        cblas_zsyrk(CblasColMajor, CblasLower, CblasNoTrans,
                    n-k-1, 1,
                    CBLAS_SADDR( alpha ), Amk, lda,
                    CBLAS_SADDR( zone  ), Akk, lda);

        /* TODO: replace by SYR but [cz]syr exists only in LAPACK */
        /* LAPACKE_zsyr_work(LAPACK_COL_MAJOR, 'l', */
        /*                   n-k-1, alpha, */
        /*                   Amk, 1, */
        /*                   Akk, lda); */

        /* Move to next Amk */
        Amk = Akk+1;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp - Computes the block static pivoting factorization of the
 * hermitian matrix n-by-n A such that A = L * D * L^t.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with Cholesky factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[in,out] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
#define MAXSIZEOFBLOCKS 64

static void core_zsytrfsp(pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        lda,
                          pastix_int_t       *nbpivot,
                          double              criteria,
                          pastix_complex64_t *work)
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, col;
    pastix_complex64_t *tmp,*tmp1,*tmp2;
    pastix_complex64_t alpha;

    /* diagonal supernode is divided into MAXSIZEOFBLOCK-by-MAXSIZEOFBLOCKS blocks */
    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    for (k=0; k<blocknbr; k++) {

        blocksize = pastix_imin(MAXSIZEOFBLOCKS, n-k*MAXSIZEOFBLOCKS);
        tmp  = A+(k*MAXSIZEOFBLOCKS)*(lda+1); /* Lk,k     */
        tmp1 = tmp  + blocksize;              /* Lk+1,k   */
        tmp2 = tmp1 + blocksize * lda;        /* Lk+1,k+1 */

        /* Factorize the diagonal block Akk*/
        core_zsytf2sp(blocksize, tmp, lda, nbpivot, criteria);

        if ((k*MAXSIZEOFBLOCKS+blocksize) < n) {

            matrixsize = n-(k*MAXSIZEOFBLOCKS+blocksize);

            /* Compute the column L(k+1:n,k) = (L(k,k)D(k,k))^{-1}A(k+1:n,k)    */
            /* 1) Compute A(k+1:n,k) = A(k+1:n,k)L(k,k)^{-T} = D(k,k)L(k+1:n,k) */
                        /* input: L(k,k) in tmp, A(k+1:n,k) in tmp1   */
                        /* output: A(k+1:n,k) in tmp1                 */
            cblas_ztrsm(CblasColMajor,
                        CblasRight, CblasLower,
                        CblasTrans, CblasUnit,
                        matrixsize, blocksize,
                        CBLAS_SADDR(zone), tmp,  lda,
                                           tmp1, lda);

            /* Compute L(k+1:n,k) = A(k+1:n,k)D(k,k)^{-1}     */
            for(col = 0; col < blocksize; col++) {
                /* copy L(k+1+col:n,k+col)*D(k+col,k+col) into work(:,col) */
                cblas_zcopy(matrixsize, tmp1+col*lda,     1,
                                        work+col*matrixsize, 1);

                                /* compute L(k+1+col:n,k+col) = A(k+1+col:n,k+col)D(k+col,k+col)^{-1} */
                alpha = 1. / *(tmp + col*(lda+1));
                cblas_zscal(matrixsize, CBLAS_SADDR(alpha),
                            tmp1+col*lda, 1);
            }

            /* Update A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - (L(k+1:n,k)*D(k,k))*L(k+1:n,k)^T */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), work, matrixsize,
                                            tmp1, lda,
                        CBLAS_SADDR(zone),  tmp2, lda);
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp1d_sytrf - Computes the LDL^H factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 * @param[in,out] work
 *          Temporary buffer used in core_zsytrfsp().
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
int core_zsytrfsp1d_sytrf( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work )
{
    pastix_int_t  ncols, stride;
    pastix_int_t  nbpivot = 0;

    ncols   = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;

    /* check if diagonal column block */
    assert( cblk->fcolnum == cblk->fblokptr->frownum );
    assert( cblk->lcolnum == cblk->fblokptr->lrownum );

    /* Factorize diagonal block (two terms version with workspace) */
    core_zsytrfsp(ncols, L, stride, &nbpivot, criteria, work);

    return nbpivot;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp1d_trsm - Apply all the trsm updates to one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
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
    SolverBlok   *blok;
    pastix_int_t  ncols, stride;
    SolverBlok   *lblk;

    ncols   = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    blok = cblk->fblokptr;   /* this diagonal block */
    lblk = cblk[1].fblokptr; /* the next diagonal block */

    /* if there are off-diagonal supernodes in the column */
    if ( blok+1 < lblk )
    {
        pastix_complex64_t *fL;
        pastix_int_t nrows, k;

        /* vertical dimension */
        nrows = stride - ncols;

        /* the first off-diagonal block in column block address */
        fL = L + blok[1].coefind;

        /* Three terms version, no need to keep L and L*D */
        cblas_ztrsm(CblasColMajor,
                    CblasRight, CblasLower,
                    CblasTrans, CblasUnit,
                    nrows, ncols,
                    CBLAS_SADDR(zone), L,  stride,
                                       fL, stride);

        for (k=0; k<ncols; k++)
        {
            pastix_complex64_t alpha;
            alpha = 1. / L[k+k*stride];
            cblas_zscal(nrows, CBLAS_SADDR(alpha), &(fL[k*stride]), 1);
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp1d_gemm - Computes the Cholesky factorization of one panel and
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
 * @param[in,out] L
 *          The pointer to the matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] work1
 *          Temporary buffer used in core_zgemdm().
 *
 * @param[in,out] work2
 *          Temporary buffer used in core_zgemdm().
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
void core_zsytrfsp1d_gemm( SolverCblk         *cblk,
                           SolverBlok         *blok,
                           SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2 )
{
    SolverBlok *iterblok;
    SolverBlok *fblok; /* facing blok */
    SolverBlok *lblok; /* first blok of panel cblk+1 */

    pastix_complex64_t *Aik, *Aij;
    pastix_int_t stride, stridefc, indblok;
    pastix_int_t dimi, dimj, dima, dimb;
    pastix_int_t ldw;
    pastix_int_t ret;
    (void)ret;

    stride = cblk->stride;
    dima = cblk_colnbr( cblk );

    /* First blok */
    indblok = blok->coefind;

    dimj = blok_rownbr( blok );
    dimi = stride - indblok;

    /* Matrix A = Aik */
    Aik = L + indblok;

    /* Compute ldw which should never be larger than SOLVE_COEFMAX */
    ldw = (dimi+1) * dima;

    /* Compute the contribution */
    ret = core_zgemdm( CblasNoTrans, CblasTrans,
                       dimi, dimj, dima,
                       1.,  Aik,   stride,
                            Aik,   stride,
                       0.,  work1, dimi,
                            L,     stride+1,
                            work2, ldw );
    assert(ret == PASTIX_SUCCESS);

    /*
     * Add contribution to facing cblk
     * A(i,i+1:n) += work1
     */

    /* Get the first block of the distant panel */
    fblok = fcblk->fblokptr;

    /* Move the pointer to the top of the right column */
    stridefc = fcblk->stride;
    C = C + (blok->frownum - fcblk->fcolnum) * stridefc;

    lblok = cblk[1].fblokptr;

    /* for all following blocks in block column */
    for (iterblok=blok; iterblok<lblok; iterblok++) {

        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;
        dimb = iterblok->lrownum - iterblok->frownum + 1;

        pastix_cblk_lock( fcblk );
        core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                     work1, dimi,
                     Aij,   stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        work1 += dimb;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp1d - Computes the Cholesky factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
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
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
int core_zsytrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           double              criteria,
                           pastix_complex64_t *work)
{
    pastix_int_t  nbpivot = core_zsytrfsp1d_sytrf(cblk, L, criteria, work);
    core_zsytrfsp1d_trsm(cblk, L);

    return nbpivot;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zsytrfsp1d - Computes the Cholesky factorization of one panel, apply all
 * the trsm updates to this panel, and apply all updates to the right with no
 * lock.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
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
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
int
core_zsytrfsp1d( SolverMatrix *solvmtx,
                 SolverCblk   *cblk,
                 double        criteria )
{
    pastix_complex64_t *L = cblk->lcoeftab;
    pastix_complex64_t *work1 = NULL;
    pastix_complex64_t *work2 = NULL;
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivot;
    pastix_int_t maxarea;

    blok = cblk->fblokptr; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    maxarea = blok_rownbr( blok ) * blok_rownbr( blok );
    blok++;
    if ( blok < lblk ) {
       for( ; blok < lblk; blok++ )
        {
            maxarea = pastix_imax( maxarea, (blok_rownbr( blok )+1) * cblk->stride );
        }
    }
    MALLOC_INTERN( work1, maxarea, pastix_complex64_t );
    MALLOC_INTERN( work2, cblk->stride * cblk_colnbr(cblk), pastix_complex64_t );

    /* if there are off-diagonal supernodes in the column */
    nbpivot = core_zsytrfsp1d_sytrf(cblk, L, criteria, work1);
    core_zsytrfsp1d_trsm(cblk, L);

    blok = cblk->fblokptr+1;
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->cblknum);

        core_zsytrfsp1d_gemm( cblk, blok, fcblk,
                              L, fcblk->lcoeftab,
                              work1, work2 );
    }

    memFree_null( work1 );
    memFree_null( work2 );
    return nbpivot;
}

