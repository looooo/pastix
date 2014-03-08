/**
 *
 * @file core_zgetrfsp.c
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
#include <math.h>
#include <cblas.h>
#include "../sopalin/sopalin_acces.h"
#include "pastix_zcores.h"

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetf2sp - Computes the sequential static pivoting LU factorization of
 * the matrix m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size stride -by- n.
 *
 * @param[in] stride
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
static void core_zgetf2sp(pastix_int_t  m,
                          pastix_int_t  n,
                          pastix_complex64_t * A,
                          pastix_int_t  stride,
                          pastix_int_t *nbpivot,
                          double criteria )
{
    pastix_int_t k, minMN;
    pastix_complex64_t *Akk, *Aik, alpha;

    minMN = pastix_imin( m, n );

    Akk = A;
    for (k=0; k<minMN; k++) {
        Aik = Akk + 1;

        if ( cabs(*Akk) < criteria ) {
            (*Akk) = (pastix_complex64_t)criteria;
            (*nbpivot)++;
        }

        /* A_ik = A_ik / A_kk, i = k+1 .. n */
        alpha = 1. / (*Akk);
        cblas_zscal(m-k-1, CBLAS_SADDR( alpha ), Aik, 1 );

        if ( k+1 < minMN ) {

            /* A_ij = A_ij - A_ik * A_kj, i,j = k+1..n */
            cblas_zgeru(CblasColMajor, m-k-1, n-k-1,
                        CBLAS_SADDR(mzone),
                        Aik,        1,
                        Akk+stride, stride,
                        Aik+stride, stride);
        }

        Akk += stride+1;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp - Computes the block static pivoting LU factorization of
 * the matrix m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size stride -by- n.
 *
 * @param[in] stride
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

static void core_zgetrfsp(pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        stride,
                          pastix_int_t       *nbpivot,
                          double              criteria)
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, tempm;
    pastix_complex64_t *Akk, *Lik, *Ukj, *Aij;

    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    Akk = A; /* Lk,k     */

    for (k=0; k<blocknbr; k++) {

        tempm = n - k * MAXSIZEOFBLOCKS;
        blocksize = pastix_imin(MAXSIZEOFBLOCKS, tempm);
        Lik = Akk + blocksize;
        Ukj = Akk + blocksize*stride;
        Aij = Ukj + blocksize;

        /* Factorize the diagonal block Akk*/
        core_zgetf2sp( tempm, blocksize, Akk, stride, nbpivot, criteria );

        matrixsize = tempm - blocksize;
        if ( matrixsize > 0 ) {

            /* Compute the column Ukk+1 */
            cblas_ztrsm(CblasColMajor,
                        CblasLeft, CblasLower,
                        CblasNoTrans, CblasUnit,
                        blocksize, matrixsize,
                        CBLAS_SADDR(zone), Akk, stride,
                                           Ukj, stride);

            /* Update Ak+1,k+1 = Ak+1,k+1 - Lk+1,k*Uk,k+1 */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), Lik, stride,
                                            Ukj, stride,
                        CBLAS_SADDR(zone),  Aij, stride);
        }

        Akk += blocksize * (stride+1);
    }
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d - Computes the LU factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize whitch Cholesky factorization. The matrix
 *          is of size stride -by- n.
 *
 * @param[in] stride
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
void core_zgetrfsp1d( pastix_complex64_t *L,
                      pastix_complex64_t *U,
                      SolverMatrix *datacode,
                      pastix_int_t c,
                      double criteria)
{
    SolverCblk *cblk;
    SolverBlok *blok;

    pastix_complex64_t *fL, *fU;
    pastix_int_t dima, dimb, stride;
    pastix_int_t fblknum, lblknum;
    pastix_int_t nbpivot = 0; /* TODO: return to higher level */
    pastix_int_t gcblk2list = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB( c ));

    cblk    = &(datacode->cblktab[c]);
    dima    = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    fblknum = cblk->bloknum;   /* block number of this diagonal block */
    lblknum = cblk[1].bloknum; /* block number of the next diagonal block */

    blok = &(datacode->bloktab[fblknum]);

    /* check if diagonal column block */
    assert( cblk->fcolnum == blok->frownum );

    /* If the block is not a leaf, we need to combine contributions on U and L */
    if ( gcblk2list != -1 )
    {
        pastix_int_t browk  = UPDOWN_LISTPTR( gcblk2list );
        pastix_int_t browk1 = UPDOWN_LISTPTR( gcblk2list+1 );

        if ( browk != browk1 ) {
            core_zaxpyt( dima, dima, 1.0,
                         U, stride,
                         L, stride );
        }
    }

    /* Factorize diagonal block (two terms version with workspace) */
    core_zgetrfsp(dima, L, stride, &nbpivot, criteria);

    /* Transpose Akk in ucoeftab */
    core_zgetro(dima, dima, L, stride, U, stride);

    /* vertical dimension */
    dimb = stride - dima;

    /* if there is an extra-diagonal bloc in column block */
    if ( fblknum+1 < lblknum )
    {
        /* first extra-diagonal bloc in column block address */
        fL = L + blok[1].coefind;
        fU = U + blok[1].coefind;

        cblas_ztrsm(CblasColMajor,
                    CblasRight, CblasUpper,
                    CblasNoTrans, CblasNonUnit,
                    dimb, dima,
                    CBLAS_SADDR(zone), L,  stride,
                                       fL, stride);

        cblas_ztrsm(CblasColMajor,
                    CblasRight, CblasUpper,
                    CblasNoTrans, CblasUnit,
                    dimb, dima,
                    CBLAS_SADDR(zone), U,  stride,
                                       fU, stride);
    }
}


void core_zgetrfsp1d_gemm(pastix_int_t cblknum,
                          pastix_int_t bloknum,
                          pastix_int_t fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U,
                          pastix_complex64_t *Cl,
                          pastix_complex64_t *Cu,
                          pastix_complex64_t *work,
                          SolverMatrix *datacode)
{
    SolverCblk *cblk  = &(datacode->cblktab[cblknum]);
    SolverCblk *fcblk = &(datacode->cblktab[fcblknum]);
    SolverBlok *blok;
    SolverBlok *fblok;

    pastix_complex64_t *Aik, *Akj, *Aij, *C;
    pastix_complex64_t *wtmp;
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
    Akj = U + indblok;

    /*
     * Compute update on L
     */
    wtmp = work;
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                 dimi, dimj, dima,
                 CBLAS_SADDR(zone),  Aik,  stride,
                                     Akj,  stride,
                 CBLAS_SADDR(zzero), wtmp, dimi  );

    /*
     * Add contribution to facing cblk
     */

    /* Get the first block of the distant panel */
    b     = fcblk->bloknum;
    fblok = &(datacode->bloktab[ b ]);

    /* Move the pointer to the top of the right column */
    stridefc = fcblk->stride;
    C = Cl + (blok->frownum - fcblk->fcolnum) * stridefc;

    lblknum = cblk[1].bloknum;

    /* for all following blocks in block column */
    for (j=bloknum; j<lblknum; j++,blok++) {

        /* Find facing bloknum */
#ifdef NAPA_SOPALIN /* ILU(k) */
        while (!(((blok->frownum >= fblok->frownum) &&
                  (blok->lrownum <= fblok->lrownum)) ||
                 ((blok->frownum <= fblok->frownum) &&
                  (blok->lrownum >= fblok->lrownum)) ||
                 ((blok->frownum <= fblok->frownum) &&
                  (blok->lrownum >= fblok->frownum)) ||
                 ((blok->frownum <= fblok->lrownum) &&
                  (blok->lrownum >= fblok->lrownum))))
#else
        while (!((blok->frownum >= fblok->frownum) &&
                 (blok->lrownum <= fblok->lrownum)))
#endif
            {
                b++; fblok++;
                assert( b < fcblk[1].bloknum );
            }

        Aij = C + fblok->coefind + blok->frownum - fblok->frownum;
        dimb = blok->lrownum - blok->frownum + 1;

        pastix_cblk_lock( fcblk );
        core_zgeadd( dimb, dimj, -1.0,
		     wtmp, dimi,
		     Aij,  stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        wtmp += dimb;
    }


    /*
     * Compute update on U
     */
    Aik = U + indblok;
    Akj = L + indblok;

    /* Restore data */
    wtmp  = work;
    blok  = &(datacode->bloktab[bloknum]);
    b     = fcblk->bloknum;
    fblok = &(datacode->bloktab[ b ]);

    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                 dimi, dimj, dima,
                 CBLAS_SADDR(zone),  Aik,  stride,
                                     Akj,  stride,
                 CBLAS_SADDR(zzero), wtmp, dimi  );

    wtmp += blok->lrownum - blok->frownum + 1;

    /*
     * Add contribution to facing cblk
     */
    C  = Cl + (blok->frownum - fcblk->fcolnum);
    Cu = Cu + (blok->frownum - fcblk->fcolnum) * stridefc;

    /*
     * Update on L part (blocks facing diagonal block)
     */
    blok++;
    for (j=bloknum+1; j<lblknum; j++, blok++) {

        /* Find facing bloknum */
#ifdef NAPA_SOPALIN /* ILU(k) */
        if (!(((blok->frownum >= fblok->frownum) &&
               (blok->lrownum <= fblok->lrownum)) ||
              ((blok->frownum <= fblok->frownum) &&
               (blok->lrownum >= fblok->lrownum)) ||
              ((blok->frownum <= fblok->frownum) &&
               (blok->lrownum >= fblok->frownum)) ||
              ((blok->frownum <= fblok->lrownum) &&
               (blok->lrownum >= fblok->lrownum))))
#else
         if (!((blok->frownum >= fblok->frownum) &&
               (blok->lrownum <= fblok->lrownum)))
#endif
             break;

        Aij = C + (blok->frownum - fblok->frownum)*stridefc;
        dimb = blok->lrownum - blok->frownum + 1;

        pastix_cblk_lock( fcblk );
        core_zaxpyt( dimj, dimb, -1.0,
                     wtmp, dimi,
                     Aij,  stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        wtmp += dimb;
    }

    C = Cu;

    /*
     * For all following blocks in block column,
     * Update is done on U directly
     */
    /* Keep updating on U */
    for (; j<lblknum; j++,blok++) {

        /* Find facing bloknum */
#ifdef NAPA_SOPALIN /* ILU(k) */
        while (!(((blok->frownum >= fblok->frownum) &&
                  (blok->lrownum <= fblok->lrownum)) ||
                 ((blok->frownum <= fblok->frownum) &&
                  (blok->lrownum >= fblok->lrownum)) ||
                 ((blok->frownum <= fblok->frownum) &&
                  (blok->lrownum >= fblok->frownum)) ||
                 ((blok->frownum <= fblok->lrownum) &&
                  (blok->lrownum >= fblok->lrownum))))
#else
        while (!((blok->frownum >= fblok->frownum) &&
                 (blok->lrownum <= fblok->lrownum)))
#endif
            {
                b++; fblok++;
                assert( b < fcblk[1].bloknum );
            }


        dimb = blok->lrownum - blok->frownum + 1;
        Aij = C + fblok->coefind + blok->frownum - fblok->frownum;

        pastix_cblk_lock( fcblk );
        core_zgeadd( dimb, dimj, -1.0,
                     wtmp, dimi,
                     Aij,  stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        wtmp += dimb;
    }
}
