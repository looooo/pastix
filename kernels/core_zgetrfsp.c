/**
 *
 * @file core_zgetrfsp.c
 *
 * PaStiX kernel routines for LU factorization.
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
#include "eztrace_module/kernels_ev_codes.h"

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
 * @brief Compute the sequential static pivoting LU factorization of the matrix
 * m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LU factorization. The matrix
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
static inline void
core_zgetf2sp( pastix_int_t        m,
               pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivot,
               double              criteria )
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
                        Akk+lda, lda,
                        Aik+lda, lda);
        }

        Akk += lda+1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the block static pivoting LU factorization of the matrix
 * m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[inout] A
 *          The matrix A to factorize with LU factorization. The matrix
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
void
core_zgetrfsp( pastix_int_t        n,
               pastix_complex64_t *A,
               pastix_int_t        lda,
               pastix_int_t       *nbpivot,
               double              criteria )
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, tempm;
    pastix_complex64_t *Akk, *Lik, *Ukj, *Aij;

    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    Akk = A; /* Lk,k     */

    for (k=0; k<blocknbr; k++) {

        tempm = n - k * MAXSIZEOFBLOCKS;
        blocksize = pastix_imin(MAXSIZEOFBLOCKS, tempm);
        Lik = Akk + blocksize;
        Ukj = Akk + blocksize*lda;
        Aij = Ukj + blocksize;

        /* Factorize the diagonal block Akk*/
        core_zgetf2sp( tempm, blocksize, Akk, lda, nbpivot, criteria );

        matrixsize = tempm - blocksize;
        if ( matrixsize > 0 ) {

            /* Compute the column Ukk+1 */
            cblas_ztrsm(CblasColMajor,
                        CblasLeft, CblasLower,
                        CblasNoTrans, CblasUnit,
                        blocksize, matrixsize,
                        CBLAS_SADDR(zone), Akk, lda,
                                           Ukj, lda);

            /* Update Ak+1,k+1 = Ak+1,k+1 - Lk+1,k*Uk,k+1 */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), Lik, lda,
                                            Ukj, lda,
                        CBLAS_SADDR(zone),  Aij, lda);
        }

        Akk += blocksize * (lda+1);
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the LU factorization of the diagonal block in a panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] U
 *          The pointer to the upper matrix storing the coefficients of the
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
cpucblk_zgetrfsp1d_getrf( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U,
                          double              criteria )
{
    pastix_int_t ncols, stride;
    pastix_int_t nbpivot = 0;

    ncols  = cblk->lcolnum - cblk->fcolnum + 1;
    stride = (cblk->cblktype & CBLK_LAYOUT_2D) ? ncols : cblk->stride;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        assert( cblk->fblokptr->LRblock[0].rk == -1 &&
                cblk->fblokptr->LRblock[1].rk == -1 );
        L = cblk->fblokptr->LRblock[0].u;
        U = cblk->fblokptr->LRblock[1].u;
        stride = ncols;

        assert( stride == cblk->fblokptr->LRblock[0].rkmax );
        assert( stride == cblk->fblokptr->LRblock[1].rkmax );
    }

    core_zgeadd( PastixTrans, ncols, ncols,
                 1.0, U, stride,
                 1.0, L, stride );

    /* Factorize diagonal block */
    start_trace_kernel( GETRF );
    core_zgetrfsp(ncols, L, stride, &nbpivot, criteria);
    stop_trace_kernel( FLOPS_ZGETRF( ncols, ncols ) );

    /* Transpose Akk in ucoeftab */
    core_zgetro(ncols, ncols, L, stride, U, stride);

    return nbpivot;
}

/**
 *******************************************************************************
 *
 * @brief Compute the LU factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[inout] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
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
 * @return The number of static pivoting performed during the diagonal block
 *         factorization.
 *
 *******************************************************************************/
int
cpucblk_zgetrfsp1d_panel( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U,
                          double              criteria,
                          const pastix_lr_t  *lowrank )
{
    pastix_int_t nbpivot;
    nbpivot = cpucblk_zgetrfsp1d_getrf(cblk, L, U, criteria);

    /*
     * We exploit the fact that the upper triangle is stored at the top of the L
     * column, and by transposition the L part of the diagonal block is
     * similarly stored in the U panel
     */
    cpucblk_ztrsmsp( PastixLCoef, PastixRight, PastixUpper,
                     PastixNoTrans, PastixNonUnit,
                     cblk, L, L, lowrank );
    cpucblk_ztrsmsp( PastixUCoef, PastixRight, PastixUpper,
                     PastixNoTrans, PastixUnit,
                     cblk, U, U, lowrank );
    return nbpivot;
}

/**
 *******************************************************************************
 *
 * @brief Perform the LU factorization of a given panel and apply all its
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
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 *******************************************************************************
 *
 * @return The number of static pivoting during factorization of the diagonal
 * block.
 *
 *******************************************************************************/
int
cpucblk_zgetrfsp1d( SolverMatrix       *solvmtx,
                    SolverCblk         *cblk,
                    double              criteria,
                    pastix_complex64_t *work )
{
    pastix_complex64_t *L = cblk->lcoeftab;
    pastix_complex64_t *U = cblk->ucoeftab;
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivot;

    nbpivot = cpucblk_zgetrfsp1d_panel(cblk, L, U, criteria, &solvmtx->lowrank);

    blok = cblk->fblokptr + 1; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    /* if there are off-diagonal supernodes in the column */
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        /* Update on L */
        cpucblk_zgemmsp( PastixLCoef, PastixUCoef, PastixTrans,
                         cblk, blok, fcblk,
                         L, U, fcblk->lcoeftab,
                         work, &solvmtx->lowrank );

        /* Update on U */
        if ( blok+1 < lblk ) {
            cpucblk_zgemmsp( PastixUCoef, PastixLCoef, PastixTrans,
                             cblk, blok, fcblk,
                             U, L, fcblk->ucoeftab,
                             work, &solvmtx->lowrank );
        }
        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }

    return nbpivot;
}
