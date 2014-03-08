/**
 *
 * @file core_zpotrfsp.c
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
#include <assert.h>

#include "common.h"
#include "pastix_zcores.h"
#include <cblas.h>

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zpotf2sp - Computes the sequential static pivoting Cholesky
 * factorization of the matrix n-by-n A = L * L^t .
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
static void core_zpotf2sp(pastix_int_t        n,
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

        /* Hermitian matrices, so imaginary part should be 0 */
        if ( creal(*Akk) < 0. )
        {
            errorPrint("Negative diagonal term\n");
            EXIT(MOD_SOPALIN, INTERNAL_ERR);
        }

        *Akk = (pastix_complex64_t)csqrt(*Akk);
        alpha = 1. / (*Akk);

        /* Scale the diagonal to compute L((k+1):n,k) */
        cblas_zscal(n-k-1, CBLAS_SADDR( alpha ), Amk, 1 );

        /* Move to next Akk */
        Akk += (lda+1);

        cblas_zher(CblasColMajor, CblasLower,
                   n-k-1, -1.,
                   Amk, 1,
                   Akk, lda);

        /* Move to next Amk */
        Amk = Akk+1;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zpotrfsp - Computes the block static pivoting Cholesky
 * factorization of the matrix n-by-n A = L * L^t .
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

static void core_zpotrfsp(pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        lda,
                          pastix_int_t       *nbpivot,
                          double              criteria)
{
    pastix_int_t k, blocknbr, blocksize, matrixsize;
    pastix_complex64_t *tmp,*tmp1,*tmp2;

    /* diagonal supernode is divided into MAXSIZEOFBLOCK-by-MAXSIZEOFBLOCKS blocks */
    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    for (k=0; k<blocknbr; k++) {

        blocksize = pastix_imin(MAXSIZEOFBLOCKS, n-k*MAXSIZEOFBLOCKS);
        tmp  = A+(k*MAXSIZEOFBLOCKS)*(lda+1); /* Lk,k     */
        tmp1 = tmp  + blocksize;              /* Lk+1,k   */
        tmp2 = tmp1 + blocksize * lda;        /* Lk+1,k+1 */

        /* Factorize the diagonal block Akk*/
        core_zpotf2sp(blocksize, tmp, lda, nbpivot, criteria);

        if ((k*MAXSIZEOFBLOCKS+blocksize) < n) {

            matrixsize = n-(k*MAXSIZEOFBLOCKS+blocksize);

            /* Compute the column L(k+1:n,k) = (L(k,k)D(k,k))^{-1}A(k+1:n,k)    */
            /* 1) Compute A(k+1:n,k) = A(k+1:n,k)L(k,k)^{-T} = D(k,k)L(k+1:n,k) */
                        /* input: L(k,k) in tmp, A(k+1:n,k) in tmp1   */
                        /* output: A(k+1:n,k) in tmp1                 */
            cblas_ztrsm(CblasColMajor,
                        CblasRight, CblasLower,
                        CblasConjTrans, CblasNonUnit,
                        matrixsize, blocksize,
                        CBLAS_SADDR(zone), tmp,  lda,
                                           tmp1, lda);

            /* Update Ak+1k+1 = Ak+1k+1 - Lk+1k * Lk+1kT */
            cblas_zherk(CblasColMajor, CblasLower, CblasNoTrans,
                        matrixsize, blocksize,
                        (double)mzone, tmp1, lda,
                        (double)zone,  tmp2, lda);
        }
    }
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zpotrfsp1d - Computes the Cholesky factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] datacode
 *          TODO
 *
 * @param[in] cblknum
 *          Index of the panel to factorize in the cblktab array.
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
int core_zpotrfsp1d( SolverMatrix       *datacode,
                     pastix_int_t        cblknum,
                     pastix_complex64_t *L,
                     double              criteria )
{
    SolverCblk   *cblk;
    SolverBlok   *blok;
    pastix_int_t  ncols, stride;
    pastix_int_t  fblknum, lblknum;
    pastix_int_t  nbpivot = 0;

    cblk    = &(datacode->cblktab[cblknum]);
    ncols   = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    fblknum = cblk->bloknum;   /* block number of this diagonal block */
    lblknum = cblk[1].bloknum; /* block number of the next diagonal block */

    /* check if diagonal column block */
    blok = &(datacode->bloktab[fblknum]);
    assert( cblk->fcolnum == blok->frownum );
    assert( cblk->lcolnum == blok->lrownum );

    /* Factorize diagonal block */
    core_zpotrfsp(ncols, L, stride, &nbpivot, criteria );

    /* if there are off-diagonal supernodes in the column */
    if ( fblknum+1 < lblknum )
    {
        pastix_complex64_t *fL;
        pastix_int_t nrows;

        /* vertical dimension */
        nrows = stride - ncols;

        /* the first off-diagonal block in column block address */
        fL = L + blok[1].coefind;

        /* Three terms version, no need to keep L and L*D */
        cblas_ztrsm(CblasColMajor,
                    CblasRight, CblasLower,
                    CblasConjTrans, CblasNonUnit,
                    nrows, ncols,
                    CBLAS_SADDR(zone), L,  stride,
                                       fL, stride);
    }

    return nbpivot;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zpotrfsp1d - Computes the Cholesky factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] datacode
 *          TODO
 *
 * @param[in] cblk
 *          The pointer to the data structure that describes the panel to be
 *          factorized. Must be at least of size 2 in order to get the first
 *          block number of the following panel (cblk[1]).
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
void core_zpotrfsp1d_gemm(SolverMatrix       *datacode,
                          pastix_int_t        cblknum,
                          pastix_int_t        bloknum,
                          pastix_int_t        fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work)
{
    SolverCblk *cblk  = &(datacode->cblktab[cblknum]);
    SolverCblk *fcblk = &(datacode->cblktab[fcblknum]);
    SolverBlok *blok;
    SolverBlok *fblok;

    pastix_complex64_t *Aik, *Aij;
    pastix_int_t lblknum;
    pastix_int_t stride, stridefc, indblok;
    pastix_int_t b, j;
    pastix_int_t dimi, dimj, dima, dimb;

    stride  = cblk->stride;
    dima = cblk->lcolnum - cblk->fcolnum + 1;


    /* First blok */
    j = bloknum;
    blok = &(datacode->bloktab[bloknum]);
    indblok = blok->coefind;

    dimj = blok->lrownum - blok->frownum + 1;
    dimi = stride - indblok;

    /* Matrix A = Aik */
    Aik = L + indblok;

    /* Compute the contribution */
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                 dimi, dimj, dima,
                 CBLAS_SADDR(zone),   Aik,  stride,
                                      Aik,  stride,
                 CBLAS_SADDR(zzero),  work, dimi );

    /*
     * Add contribution to facing cblk
     * A(i,i+1:n) += work1
     */

    /* Get the first block of the distant panel */
    b     = fcblk->bloknum;
    fblok = &(datacode->bloktab[ b ]);

    /* Move the pointer to the top of the right column */
    stridefc = fcblk->stride;
    C = C + (blok->frownum - fcblk->fcolnum) * stridefc;

    lblknum = cblk[1].bloknum;

    /* for all following blocks in block column */
    for (j=bloknum; j<lblknum; j++,blok++) {

        /* Find facing bloknum */
        while (!is_block_inside_fblock( blok, fblok ))
        {
            b++; fblok++;
            assert( b < fcblk[1].bloknum );
        }

        Aij = C + fblok->coefind + blok->frownum - fblok->frownum;
        dimb = blok->lrownum - blok->frownum + 1;

        pastix_cblk_lock( fcblk );
        core_zgeadd( dimb, dimj, -1.0,
                     work, dimi,
                     Aij,  stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        work += dimb;
    }
}
