/**
 *
 * @file core_zgeadd.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#include "common.h"
#include "pastix_zcores.h"
#include <cblas.h>

/**
 ******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 *  core_zgeadd adds to matrices together.
 *
 *       B <- alpha * op(A)  + B
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg CblasNoTrans   :  No transpose, op( A ) = A;
 *         @arg CblasTrans     :  Transpose, op( A ) = A';
 *         @arg CblasConjTrans :  Conjugate Transpose, op( A ) = conj(A').
 *
 * @param[in] M
 *          Number of rows of the matrix B.
 *          Number of rows of the matrix A, if trans == CblasNoTrans, number of
 *          columns of A otherwise.
 *
 * @param[in] N
 *          Number of columns of the matrix B.
 *          Number of columns of the matrix A, if trans == CblasNoTrans, number
 *          of rows of A otherwise.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans == CblasNoTrans, LDA-by-M,
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,K).
 *          K = M if trans == CblasNoTrans, K = N otherwise.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval 1, not yet implemented
 *
 ******************************************************************************/
int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB)
{
    int j;

#if defined(PASTIX_DEBUG)
    if (M < 0) {
        //coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        //coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < pastix_imax(1,M)) && (M > 0) ) {
        //coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < pastix_imax(1,M)) && (M > 0) ) {
        //coreblas_error(7, "Illegal value of LDB");
        return -7;
    }
#endif

    if (trans == CblasNoTrans) {
        if (M == LDA && M == LDB)
            cblas_zaxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
        else {
            for (j = 0; j < N; j++)
                cblas_zaxpy(M, CBLAS_SADDR(alpha), A + j*LDA, 1, B + j*LDB, 1);
        }
    }
    else if (trans == CblasTrans ) {
        for (j = 0; j < N; j++)
            cblas_zaxpy(M, CBLAS_SADDR(alpha), A + j, LDA, B + j*LDB, 1);

    }
    /* trans == CblasConjTrans */
    else {
        return -10;
    }

    return PASTIX_SUCCESS;
}


/**
 ******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 *  core_zgeaddsp1d adds two column blocks together.
 *
 *       cblk2 <- cblk1  + cblk2
 *
 *******************************************************************************
 *
 * @param[in] cblk1
 *          The pointer to the data structure that describes the panel to add.
 *          Next column blok must be accessible through cblk1[1].
 *
 * @param[in] cblk2
 *          The pointer to the data structure that describes the panel in which
 *          we add.
 *          Next column blok must be accessible through cblk2[1].
 *
 * @param[in] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] Cl
 *          The pointer to the lower matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width. Ignored if
 *          NULL.
 *
 *
 * @param[in,out] Cu
 *          The pointer to the upper matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
int
core_zgeaddsp1d(SolverCblk * cblk1,
                SolverCblk * cblk2,
                pastix_complex64_t * L,
                pastix_complex64_t * Cl,
                pastix_complex64_t * U,
                pastix_complex64_t * Cu) {
    SolverBlok *iterblok;
    SolverBlok *firstblok;
    SolverBlok *lastblok;
    SolverBlok *fblok;
    pastix_int_t ncol1 = cblk1->lcolnum - cblk1->fcolnum + 1;
    pastix_complex64_t *ga, *gb;

    firstblok = cblk1->fblokptr;
    lastblok  = cblk1[1].fblokptr;
    fblok = cblk2->fblokptr;

    assert(cblk1->fcolnum >= cblk2->fcolnum);

    for (iterblok = firstblok; iterblok < lastblok; iterblok++) {
        pastix_int_t nrow;
        /* Find facing bloknum */
        while (!is_block_inside_fblock( iterblok, fblok )) {
            fblok++;
            assert( fblok < cblk2[1].fblokptr );
        }
        ga = L + iterblok->coefind;
        gb = Cl + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
            fblok->coefind +
            iterblok->frownum - fblok->frownum;
        nrow = iterblok->lrownum - iterblok->frownum + 1;
        core_zgeadd( CblasNoTrans,
                     nrow, ncol1, -1.0,
		     ga, cblk1->stride,
		     gb, cblk2->stride );
        if (U != NULL) {
            ga = U + iterblok->coefind;
            gb = Cu + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
                fblok->coefind +
                iterblok->frownum - fblok->frownum;
            core_zgeadd( CblasNoTrans,
                         nrow, ncol1, -1.0,
                         ga, cblk1->stride,
                         gb, cblk2->stride );
        }
    }
    return PASTIX_SUCCESS;
}
